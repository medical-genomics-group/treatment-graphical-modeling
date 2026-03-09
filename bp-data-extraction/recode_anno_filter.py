import hail as hl
import json
import pyspark
from pyspark.sql.functions import col
import dxpy
import dxdata
import random
from bokeh.io import show, output_notebook
from bokeh.layouts import gridplot
from collections import defaultdict
import pandas as pd
import scipy.stats as stats 

def anno_bp(full):
    
    with open('recoding_dict.json') as json_file:
        recoding = json.load(json_file)

    recoding = recoding['bp']
    systems = recoding['systems']

    drug_names = hl.literal(list(recoding['drugs']['drug_name'].keys()))
    drug_name_idx = hl.map(lambda x: full.info.lower().contains(x), drug_names).index(True)
    drug_codes = hl.literal(list(recoding['drugs']['drug_name'].values()))

    full = full.annotate(
        recode = hl.if_else(
            full.system == 'read_2',
            hl.dict(systems['read_2']).get(full.code, hl.missing(hl.tstr)),
            hl.if_else(
                full.system =='bnf',
                hl.dict(systems['bnf_code']).get(full.code, hl.missing(hl.tstr)),
                hl.if_else(
                    full.system == 'dmd',
                    hl.dict(systems['dmd_code']).get(full.code, hl.missing(hl.tstr)),
                    hl.if_else(
                        full.system == 'read_3',
                        full.code,
                        hl.null(hl.tstr)
                    )
                )
            )
        ),
        recode_from_name = hl.dict(
            recoding['drugs']['drug_name']
        ).get(
            drug_names[drug_name_idx]
        )
    )

    full = full.annotate(
        recode = hl.coalesce(
            full.recode, full.recode_from_name
        )
    )

    to_text = recoding['drugs']['to_text'].copy()
    to_text.update(recoding['procedures']['to_text'])

    full = full.annotate(
        recode_text = hl.dict(to_text).get(full.recode)
    )
    
    full = full.drop(full.recode_from_name)
    
    return(full)

def filter_bp(full):
    
    import json

    with open('recoding_dict.json') as json_file:
        recoding = json.load(json_file)

    recoding = recoding['bp']
    to_text = recoding['drugs']['to_text'].copy()
    to_text.update(recoding['procedures']['to_text'])

    bp_filter = hl.literal(list(to_text.keys()))
    bp_mes = hl.literal(['O/E - blood pressure',
          'O/E-blood pressure reading NOS',
          'Lying blood pressure reading',
          'Standing blood pressure reading',
          'Sitting blood pressure reading',
          'Lying blood pressure reading',
          'ABP - Arterial blood pressure',
          'Systolic blood pressure',
          'DBP - Diastolic blood pressure',
          'Standing systolic blood pressure',
          'Standing diastolic blood pressure',
          'Sitting systolic blood pressure',
          'Sitting diastolic blood pressure',
          'Lying systolic blood pressure',
          'Lying diastolic blood pressure'])
    
    full = full.filter(
        bp_filter.contains(full.recode)
    )
    
    full = full.filter(
        hl.if_else(
            (full.source == 'GP'),
            bp_mes.contains(full.recode_text),
            hl.if_else(
                (full.source == 'PRESC'),
                True,
                False
            )
        )
    )
    
    return full

def hail_answer_q(ht,
                  agg_expr,
                  filter_expr=True):
    
    """this function will return the aggregation result as hail 
    table works best if the result is a simple dict"""
    
    res = ht.aggregate(
        hl.agg.filter(
            filter_expr,
            agg_expr
        )
    )
        
    res_ht = hl.utils.range_table(len(res))
    
    res_ht = res_ht.annotate(
        field_name = hl.literal(list(res.keys()))[res_ht.idx],
        field_value = hl.literal(list(res.values()))[res_ht.idx]
    )
    
    res_ht = res_ht.key_by(res_ht.field_name)
    res_ht = res_ht.drop(res_ht.idx)
    
    return res_ht 

def parse_bp(bp_meas):
    
    part_1 = bp_meas.info.replace(r'[^\d,\._]', '').split('_')[0]
    part_2 = bp_meas.info.replace(r'[^\d,\._]', '').split('_')[1]
    
    #note: info field had three parts, but it was manually confirmed that the last part had no measurments

    bp_meas = bp_meas.annotate(
        info_1 = part_1,
        info_2 = part_2,
    ) 

    bp_meas = bp_meas.annotate(
        info_1 = hl.if_else(
            (bp_meas.info_1 == ""),
            hl.missing(hl.tint),
            hl.int(hl.float(bp_meas.info_1))
        ),
        info_2 = hl.if_else(
            (bp_meas.info_2 == ""),
            hl.missing(hl.tint),
            hl.int(hl.float(bp_meas.info_2))
        ),
    )
    
    dias_filter = (bp_meas.recode.upper().contains('DIASTOLIC'))
                                                      
    bp_meas = bp_meas.annotate(
        systolic = hl.if_else(
            hl.is_defined(bp_meas.info_1) & hl.is_defined(bp_meas.info_2),
            hl.max(bp_meas.info_1, bp_meas.info_2),
            #hl.missing(hl.tint) # more conservative filter
            hl.if_else(
                dias_filter,
                hl.missing(hl.tint),
                hl.coalesce(bp_meas.info_1, bp_meas.info_2)
            )   
        ),
        diastolic = hl.if_else(
            hl.is_defined(bp_meas.info_1) & hl.is_defined(bp_meas.info_2),
            hl.min(bp_meas.info_1, bp_meas.info_2),
            #hl.missing(hl.tint)  # more conservative filter
            hl.if_else(
                dias_filter,
                hl.coalesce(bp_meas.info_1, bp_meas.info_2),
                hl.missing(hl.tint)
            )
        ) 
    )

   #remove impossible measurments
    bp_meas = bp_meas.annotate(
        systolic = hl.if_else(
            ((bp_meas.systolic <= 1) | (bp_meas.systolic >= 1000)),
            hl.missing(hl.tint),
            bp_meas.systolic
        ),
        diastolic = hl.if_else(
            ((bp_meas.diastolic <= 1) | (bp_meas.diastolic >= 1000)),
            hl.missing(hl.tint),
            bp_meas.diastolic
        )
    ) 
    
    systolic_stats = bp_meas.aggregate(
        hl.agg.stats(bp_meas.systolic)
    )
        
    systolic_sd =  systolic_stats['stdev']
    systolic_mean = systolic_stats['mean']
    
    systolic_lowest = systolic_stats['mean'] - (8*systolic_stats['stdev'])
    systolic_highest = systolic_stats['mean'] + (8*systolic_stats['stdev'])
    
    diastolic_stats = bp_meas.aggregate(
        hl.agg.stats(bp_meas.diastolic)
    )
    
    diastolic_sd =  diastolic_stats['stdev']
    diastolic_mean = diastolic_stats['mean']
    
    diastolic_lowest = diastolic_stats['mean'] - (8*diastolic_stats['stdev'])
    diastolic_highest = diastolic_stats['mean'] + (8*diastolic_stats['stdev'])
        
    #remove improbable measurments
    bp_meas = bp_meas.annotate(
        systolic = hl.if_else(
            ((bp_meas.systolic <= systolic_lowest) | (bp_meas.systolic >= systolic_highest)),
            hl.missing(hl.tint),
            bp_meas.systolic
        ),
        diastolic = hl.if_else(
            ((bp_meas.diastolic <= diastolic_lowest) | (bp_meas.diastolic >= diastolic_highest)),
            hl.missing(hl.tint),
            bp_meas.diastolic
        )
    ) 
    
    #filter missing measurments
    bp_meas = bp_meas.filter(
        hl.is_defined(
            hl.coalesce(bp_meas.systolic, bp_meas.diastolic)
        )
    )
    
    bp_meas = bp_meas.select(
        bp_meas.eid,
        bp_meas.source,
        bp_meas.date,
        bp_meas.systolic,
        bp_meas.diastolic   
    )
    
    return(bp_meas)


def get_bp_ases(spark):
    
    "extracts BP measurments from the assesment center to later add it to all measurments"
    
    dispensed_database_name = dxpy.find_one_data_object(
        classname="database", name="app*", folder="/", name_mode="glob", describe=True
    )["describe"]["name"]

    spark.sql("USE " + dispensed_database_name)
    
    ases = spark.sql("""
          SELECT gp.eid AS eid,
          'ASCEN' as source,
          '246..' as code,
          participant_0001.p53_i0 as date,
          'read_3' as system,
          CONCAT (COALESCE(participant_0002.p93_i0_a0, participant_0008.p4080_i0_a0),
                  '_',
                  COALESCE(participant_0002.p93_i0_a1, participant_0008.p4080_i0_a1),
                  '_'
                  ) as systolic_info,
          CONCAT (COALESCE(participant_0002.p94_i0_a0, participant_0008.p4079_i0_a0),
                  '_',
                  COALESCE(participant_0002.p94_i0_a1, participant_0008.p4079_i0_a1),
                  '_'
                  ) as diastolic_info,
          '246..' as recode,
          'blood pressure' as recode_text
          FROM gp_registrations AS gp
          JOIN participant_0001 ON gp.eid = participant_0001.eid
          JOIN participant_0002 ON gp.eid = participant_0002.eid
          JOIN participant_0008 ON gp.eid = participant_0008.eid
          
          UNION
          
          SELECT gp.eid AS eid,
          'ASCEN' as source,
          '246..' as code,
          participant_0001.p53_i1 as date,
          'read_3' as system,
          CONCAT (COALESCE(participant_0002.p93_i1_a0, participant_0008.p4080_i1_a0),
                  '_',
                  COALESCE(participant_0002.p93_i1_a1, participant_0008.p4080_i1_a1),
                  '_'
                  ) as systolic_info,
          CONCAT (COALESCE(participant_0002.p94_i1_a0, participant_0008.p4079_i1_a0),
                  '_',
                  COALESCE(participant_0002.p94_i1_a1, participant_0008.p4079_i1_a1),
                  '_'
                  ) as diastolic_info,
          '246..' as recode,
          'blood pressure' as recode_text
          
          FROM gp_registrations AS gp
          JOIN participant_0001 ON gp.eid = participant_0001.eid
          JOIN participant_0002 ON gp.eid = participant_0002.eid
          JOIN participant_0008 ON gp.eid = participant_0008.eid
          WHERE participant_0001.p53_i1 IS NOT NULL
          
          UNION
          SELECT gp.eid AS eid,
          'ASCEN' as source,
          '246..' as code,
          participant_0001.p53_i2 as date,
          'read_3' as system,
          CONCAT (COALESCE(participant_0002.p93_i2_a0, participant_0008.p4080_i2_a0),
                  '_',
                  COALESCE(participant_0002.p93_i2_a1, participant_0008.p4080_i2_a1),
                  '_'
                  ) as systolic_info,
          CONCAT (COALESCE(participant_0002.p94_i2_a0, participant_0008.p4079_i2_a0),
                  '_',
                  COALESCE(participant_0002.p94_i2_a1, participant_0008.p4079_i2_a1),
                  '_'
                  ) as diastolic_info,
          '246..' as recode,
          'blood pressure' as recode_text
                  
          FROM gp_registrations AS gp
          JOIN participant_0001 ON gp.eid = participant_0001.eid
          JOIN participant_0002 ON gp.eid = participant_0002.eid
          JOIN participant_0008 ON gp.eid = participant_0008.eid
          
          WHERE participant_0001.p53_i2 IS NOT NULL
          
          
          UNION
          SELECT gp.eid AS eid,
          'ASCEN' as source,
          '246..' as code,
          participant_0001.p53_i3 as date,
          'read_3' as system,
          CONCAT (COALESCE(participant_0002.p93_i3_a0, participant_0008.p4080_i3_a0),
                  '_',
                  COALESCE(participant_0002.p93_i3_a1, participant_0008.p4080_i3_a1),
                  '_'
                  ) as systolic_info,
          CONCAT (COALESCE(participant_0002.p94_i3_a0, participant_0008.p4079_i3_a0),
                  '_',
                  COALESCE(participant_0002.p94_i3_a1, participant_0008.p4079_i3_a1),
                  '_'
                  ) as diastolic_info,
          '246..' as recode,
          'blood pressure' as recode_text
          FROM gp_registrations AS gp
          JOIN participant_0001 ON gp.eid = participant_0001.eid
          JOIN participant_0002 ON gp.eid = participant_0002.eid
          JOIN participant_0008 ON gp.eid = participant_0008.eid
          WHERE participant_0001.p53_i3 IS NOT NULL
          """)
    
    ases = (
            ases
            .withColumn('date', col('date').cast('string'))
        )

    ases_ht = hl.Table.from_spark(ases)

    cond_1 = hl.is_defined(
        ases_ht.systolic_info.split('_')[0]
    ) & hl.is_defined(
        ases_ht.systolic_info.split('_')[1]
    )

    mean_1 = hl.mean(
        hl.set(
            [hl.int(ases_ht.systolic_info.split('_')[0]),
             hl.int(ases_ht.systolic_info.split('_')[1])]
        )
    )

    cond_2 = hl.is_defined(ases_ht.systolic_info.split('_')[0])
    cond_3 = hl.is_defined(ases_ht.systolic_info.split('_')[1])

    cond_4 = hl.is_defined(
        ases_ht.diastolic_info.split('_')[0]
    ) & hl.is_defined(
        ases_ht.diastolic_info.split('_')[1]
    )

    mean_4 = hl.mean(
        hl.set(
            [hl.int(ases_ht.diastolic_info.split('_')[0]),
             hl.int(ases_ht.diastolic_info.split('_')[1])]
        )
    )

    cond_5 = hl.is_defined(ases_ht.diastolic_info.split('_')[0])
    cond_6 = hl.is_defined(ases_ht.diastolic_info.split('_')[1])

    ases_ht = ases_ht.transmute(
        systolic =
        hl.if_else(
            cond_1,
            mean_1,
            hl.if_else(
                cond_2,
                hl.int(ases_ht.systolic_info.split('_')[0]),
                hl.if_else(
                    cond_3,
                    hl.int(ases_ht.systolic_info.split('_')[1]),
                    hl.missing(hl.tint)
                )
            )
        ),
        diastolic =
        hl.if_else(
            cond_4,
            mean_4,
            hl.if_else(
                cond_5,
                hl.int(ases_ht.diastolic_info.split('_')[0]),
                hl.if_else(
                    cond_6,
                    hl.int(ases_ht.diastolic_info.split('_')[1]),
                    hl.missing(hl.tint)
                )
            )
        )
    )
    
    ases_ht = ases_ht.annotate(
        systolic = hl.int(ases_ht.systolic),
        diastolic = hl.int(ases_ht.diastolic)
    )

    ases_ht = ases_ht.select(
        ases_ht.eid,
        ases_ht.source,
        ases_ht.date,
        ases_ht.systolic,
        ases_ht.diastolic
    )
    
    return ases_ht


def extract_age(bp_meas, spark):

    dispensed_database_name = dxpy.find_one_data_object(
        classname="database", name="app*", folder="/", name_mode="glob", describe=True
    )["describe"]["name"]

    spark.sql("USE " + dispensed_database_name)

    ages = spark.sql("""
        SELECT gp.eid AS eid,
        participant_0001.p34

        FROM gp_registrations AS gp
        JOIN participant_0001 ON gp.eid = participant_0001.eid 
              """)

    ages_ht = hl.Table.from_spark(ages)
    ages_ht = ages_ht.key_by(ages_ht.eid)

    date = hl.int(bp_meas.date.split('-')[0])
    bp_meas = bp_meas.annotate(
        age = date - ages_ht[bp_meas.eid].p34
    )
    
    return bp_meas

def filter_meas(bp, spark):
    
    bp_meas = bp.filter(bp.source == 'GP')
    bp_meas = parse_bp(bp_meas)

    ases = get_bp_ases(spark)
    bp_meas = bp_meas.union(ases)

    bp_meas = extract_age(bp_meas, spark)
    
    bp_meas = bp_meas.group_by(
        bp_meas.eid,
        bp_meas.source,
        bp_meas.date,
        bp_meas.age
    ).aggregate(
        systolic = hl.agg.filter(
            (bp_meas.systolic != hl.float('nan')),
            hl.mean(hl.agg.collect(bp_meas.systolic))
        ),
        diastolic = hl.agg.filter(
            (bp_meas.diastolic != hl.float('nan')),
            hl.mean(hl.agg.collect(bp_meas.diastolic))
        )
    )
    
    bp_meas = bp_meas.annotate(
        systolic = hl.if_else(
            (bp_meas.systolic == hl.float('nan')),
            hl.missing(hl.tfloat),
            bp_meas.systolic
        ),
        diastolic = hl.if_else(
            (bp_meas.diastolic == hl.float('nan')),
            hl.missing(hl.tfloat),
            bp_meas.diastolic
        )
    )

    bp_meas = bp_meas.key_by()
    
    return(bp_meas)

def to_days(date_field):
    year = hl.int(date_field.split('-')[0])
    month = hl.int(date_field.split('-')[1])
    day = hl.int(date_field.split('-')[2])
    result = (year*364.24)+(month*30.44)+day
    return(result)

def show_stats(bp_meas): 
    
    for_sys = bp_meas.filter(hl.is_nan(bp_meas.systolic), keep = False)
    
    bp_stats_systolic = for_sys.aggregate(
        hl.agg.group_by(
            for_sys.source,
            hl.agg.stats(for_sys.systolic)
        )
    )
    
    for_dias = bp_meas.filter(hl.is_nan(bp_meas.diastolic), keep = False)

    bp_stats_diastolic = for_dias.aggregate(
        hl.agg.group_by(
            for_dias.source,
            hl.agg.stats(for_dias.diastolic)
        )
    )

    systolic_gp_hist = bp_meas.aggregate(
        hl.agg.filter(
            bp_meas.source == 'GP',
            hl.agg.hist(bp_meas.systolic, 0, 250, 250)
        )
    )

    systolic_ac_hist = bp_meas.aggregate(
        hl.agg.filter(
            bp_meas.source == 'ASCEN',
            hl.agg.hist(bp_meas.systolic, 0, 250, 250)
        )
    )

    p = hl.plot.histogram(
        systolic_gp_hist,
        title='systolic gp hist')

    p2 = hl.plot.histogram(
        systolic_ac_hist,
        title='systolic ac hist')


    diastolic_gp_hist = bp_meas.aggregate(
        hl.agg.filter(
            bp_meas.source == 'GP',
            hl.agg.hist(bp_meas.diastolic, 0, 250, 250)
        )
    )

    diastolic_ac_hist = bp_meas.aggregate(
        hl.agg.filter(
            bp_meas.source == 'ASCEN',
            hl.agg.hist(bp_meas.diastolic, 0, 250, 250)
        )
    )

    p3 = hl.plot.histogram(
        diastolic_gp_hist,
        title='diastolic gp hist')

    p4 = hl.plot.histogram(
        diastolic_ac_hist,
        title='diastolic ac hist')


    bp_meas = bp_meas.key_by(bp_meas.eid)
    per_part = bp_meas.collect_by_key()

    per_part.aggregate(hl.agg.stats(hl.len(per_part.values)))

    per_part_hist = per_part.aggregate(
        hl.agg.hist(
            hl.len(per_part.values),
            0,
            500,
            500
        )
    ) # for most people we have ony 1 measurment

    p5 = hl.plot.histogram(
        per_part_hist,
        title='number of measurments per participant hist') 

    print("stats for systolic measurments", bp_stats_systolic)
    print("stats for diastolic measurments", bp_stats_diastolic)

    show(gridplot([p, p2, p3, p4], ncols=2))
    
    show(p5)
    
    
    
def anno_presc(bp, bp_meas, db_uri, include_day_zero):
    
    bp_meas = bp_meas.annotate(
            date = hl.int(to_days(bp_meas.date))
        )

    bp_drugs = bp.filter(bp.source == 'PRESC')
    bp_drugs = bp_drugs.key_by(bp_drugs.eid)
    bp_meas = bp_meas.key_by(bp_meas.eid)

    bp_drugs = bp_drugs.annotate(
            day = hl.int(to_days(bp_drugs.date)),
            meas_date = bp_meas.index(bp_drugs.eid, all_matches=True)['date']                
        )

    bp_drugs = bp_drugs.select(
            bp_drugs.day,
            bp_drugs.meas_date,
            bp_drugs.recode_text
        )

    bp_drugs = bp_drugs.persist()
    bp_drugs = bp_drugs.explode(bp_drugs.meas_date)

    if include_day_zero:
        bp_drugs = bp_drugs.filter(
            (bp_drugs.day >= bp_drugs.meas_date - 90) &
            (bp_drugs.day <= bp_drugs.meas_date) # 0 to 90 days for 0 to X windows
            )
    else:
        bp_drugs = bp_drugs.filter(
            (bp_drugs.day >= bp_drugs.meas_date - 90) &
            (bp_drugs.day < bp_drugs.meas_date)  # 1 to 90 days for 0 to X windows
        )
            

    bp_drugs = bp_drugs.annotate(
            date = bp_drugs.meas_date
        )

    bp_drugs = bp_drugs.annotate(
        diff_days = bp_drugs.meas_date - bp_drugs.day
    )

    diff_map_tbl = (
        bp_drugs
        .group_by(bp_drugs.eid, bp_drugs.meas_date)
        .aggregate(
            diff_map = hl.agg.group_by(
                bp_drugs.recode_text,
                hl.agg.min(bp_drugs.diff_days)
            )
        )
        .rename({'meas_date': 'date'})
        .key_by('eid', 'date')
    )

    bp_meas = bp_meas.key_by(bp_meas.eid, bp_meas.date)
    bp_meas = bp_meas.annotate(
        diff_map = diff_map_tbl.index(bp_meas.key).diff_map
    )

    unique_drugs = list(
        bp_drugs.aggregate(hl.agg.collect_as_set(bp_drugs.recode_text))
    )

    bp_meas = bp_meas.annotate(**{
        drug: hl.cond(
            hl.is_defined(bp_meas.diff_map),
            bp_meas.diff_map.get(drug),          
            hl.missing(hl.tint32)                   
        )
        for drug in unique_drugs
    })

    sources = ['GP', 'ASCEN']

    bp_meas = bp_meas.annotate(**{s:
                                  hl.int(bp_meas.source.contains(s))
                                  for s in sources})

    bp_meas = bp_meas.drop(bp_meas.diff_map, bp_meas.source)

    tb_name = "bp-for-export.ht"
    url = f"dnax://{db_uri}/{tb_name}"
    
    bp_meas.write(url, overwrite = True)
    bp_meas = hl.read_table(url)

    return(bp_meas)


def get_cvd(full, spark):
    
    with open('recoding_dict.json') as json_file:
        recoding = json.load(json_file)

    recoding = recoding['cvd']

    full = full.annotate(
        recode = hl.if_else(
            full.system == 'read_2',
            hl.dict(recoding['read_2']).get(full.code, hl.missing(hl.tstr)),
            hl.if_else(
                full.system =='read_3',
                hl.dict(recoding['read_3']).get(full.code, hl.missing(hl.tstr)),
                hl.if_else(
                    full.system =='ICD10',
                    hl.dict(recoding['ICD10']).get(full.code, hl.missing(hl.tstr)),
                    hl.if_else(
                        full.system == 'ICD9',
                        hl.dict(recoding['ICD9']).get(full.code, hl.missing(hl.tstr)),
                        hl.missing(hl.tstr)
                    )
                )
            )
        )
    )

    full = full.annotate(
        recode_text = hl.dict(recoding['to_text']).get(full.recode)
    )

    heart_disease = full.filter(
        hl.is_defined(
            full.recode
        )
    )

    heart_disease = extract_age(heart_disease, spark)
    heart_disease = heart_disease.select(heart_disease.eid, heart_disease.age)

    heart_disease = heart_disease.order_by(hl.asc('eid'), hl.asc('age'))
    
    heart_disease = (
    heart_disease.group_by(heart_disease.eid)
                .aggregate(age = hl.agg.min(heart_disease.age)))
    
    return(heart_disease)

def extract_ldl(spark):
    
    dispensed_database_name = dxpy.find_one_data_object(
        classname="database", name="app*", folder="/", name_mode="glob", describe=True
    )["describe"]["name"]
    spark.sql("USE " + dispensed_database_name)

    gp_df = spark.sql("SELECT * FROM gp_clinical")
    for col_name, data_type in gp_df.dtypes:
        if data_type == "date":
            gp_df = gp_df.withColumn(col_name, gp_df[col_name].cast("string"))

    gp_num_df = spark.sql("""
        SELECT 
            eid, 
            p42040 AS gp_num
        FROM participant_0077
    """)

    ldl_df = spark.sql("""
        SELECT 
            eid, 
            p30780_i0 AS ldl_i0,
            p30780_i1 AS ldl_i1
        FROM participant_0073
    """)

    dates_df = spark.sql("""
        SELECT 
            eid, 
            CAST(p53_i0 AS STRING) AS date_i0,
            CAST(p53_i1 AS STRING) AS date_i1
        FROM participant_0001
    """)

    ldl_ht    = hl.Table.from_spark(ldl_df)
    dates_ht  = hl.Table.from_spark(dates_df)
    gp_num_ht = hl.Table.from_spark(gp_num_df)
    gp_ht     = hl.Table.from_spark(gp_df)

    long_rows = ldl_ht.select(
        eid = ldl_ht.eid,
        ldl_data = hl.array([
            hl.struct(visit_num='0', ldl=ldl_ht.ldl_i0),
            hl.struct(visit_num='1', ldl=ldl_ht.ldl_i1),
        ])
    )
    long_table = long_rows.explode('ldl_data').select(
        eid       = long_table.eid,
        visit_num = long_table.ldl_data.visit_num,
        ldl       = long_table.ldl_data.ldl
    )
    long_table = long_table.filter(hl.is_defined(long_table.ldl))

    dates_ht = dates_ht.key_by()         
    long_table = long_table.key_by('eid')

    dates_dict_ht = dates_ht.select(
        'eid',
        dates = hl.dict([('date_i0', dates_ht.date_i0),
                         ('date_i1', dates_ht.date_i1)])
    ).key_by('eid')

    combined_table = long_table.annotate(
        _d = dates_dict_ht[long_table.eid],
        date = dates_dict_ht[long_table.eid]['dates'][('date_i' + long_table.visit_num)]
    ).drop('_d', 'visit_num')

    assestment_centre_ht = combined_table.annotate(source = "ASCEN").persist()

    gp_record_counts = gp_ht.group_by(gp_ht.eid).aggregate(record_count = hl.agg.count())
    gp_num_ht = gp_num_ht.key_by('eid')
    gp_record_counts = gp_record_counts.key_by('eid')

    comparison_table = gp_num_ht.annotate(
        counted_records = gp_record_counts[gp_num_ht.eid].record_count
    ).annotate(match = comparison_table.gp_num == comparison_table.counted_records)

    mismatches = comparison_table.filter(~comparison_table.match)
    _mismatched_count = mismatches.count()
    if _mismatched_count > 0:
        print(f"WARNING, wrong number of records in GP table: {_mismatched_count}")

    read3_ldl_read_codes = ['44P6.', '44PI.', '44PD.', 'XaIp4']
    read2_ldl_read_codes = ['44P6.', '44PI',  '44PD.']
    read3_set = hl.set(read3_ldl_read_codes)
    read2_set = hl.set(read2_ldl_read_codes)

    gp_filtered_ht = gp_ht.filter(
        read2_set.contains(gp_ht.read_2) | read3_set.contains(gp_ht.read_3)
    )

    float_pattern = r'^-?\d+(\.\d+)?$'
    gp_annotated_ht = gp_filtered_ht.annotate(
        ldl = hl.cond(
            gp_filtered_ht.value1.matches(float_pattern),
            hl.float64(gp_filtered_ht.value1),
            hl.cond(
                gp_filtered_ht.value2.matches(float_pattern),
                hl.float64(gp_filtered_ht.value2),
                hl.null(hl.tfloat64)
            )
        )
    ).drop('data_provider', 'read_2', 'read_3', 'value1', 'value2', 'value3') \
     .rename({'event_dt': 'date'}) \
     .annotate(source = "GP")

    gp_annotated_ht     = gp_annotated_ht.key_by('eid')
    assestment_centre_ht = assestment_centre_ht.key_by('eid')

    final_ldl_ht = assestment_centre_ht.union(gp_annotated_ht, unify=True).persist()

    MIN_LDL_LVL = 0.26
    MAX_LDL_LVL = 10.30
    final_ldl_ht = final_ldl_ht.filter(
        (final_ldl_ht.ldl >= MIN_LDL_LVL) & (final_ldl_ht.ldl <= MAX_LDL_LVL)
    )

    return final_ldl_ht

def anno_presc_ldl(ldl_meas, bp, db_uri):

    ldl_meas = ldl_meas.annotate(date = hl.int(to_days(ldl_meas.date)))

    bp_drugs = bp.filter(bp.source == 'PRESC')
    bp_drugs = bp_drugs.key_by(bp_drugs.eid)
    ldl_meas = ldl_meas.key_by(ldl_meas.eid)

    bp_drugs = bp_drugs.annotate(
        day = hl.int(to_days(bp_drugs.date)),
        meas_date = ldl_meas.index(bp_drugs.eid, all_matches=True)['date']
    ).select('day', 'meas_date', 'recode_text')

    bp_drugs = bp_drugs.persist()
    bp_drugs = bp_drugs.explode(bp_drugs.meas_date)

    bp_drugs = bp_drugs.filter(
        (bp_drugs.day >= bp_drugs.meas_date - 90) &
        (bp_drugs.day <= bp_drugs.meas_date)
    )

    bp_drugs = bp_drugs.annotate(
        diff_days = bp_drugs.meas_date - bp_drugs.day
    )

    diff_map_tbl = (
        bp_drugs
        .group_by(bp_drugs.eid, bp_drugs.meas_date)
        .aggregate(
            diff_map = hl.agg.group_by(
                bp_drugs.recode_text,
                hl.agg.min(bp_drugs.diff_days)
            )
        )
        .rename({'meas_date': 'date'})
        .key_by('eid', 'date')
    )

    ldl_meas = ldl_meas.key_by(ldl_meas.eid, ldl_meas.date)
    ldl_meas = ldl_meas.annotate(
        diff_map = diff_map_tbl.index(ldl_meas.key).diff_map
    )

    unique_drugs = list(
        bp_drugs.aggregate(hl.agg.collect_as_set(bp_drugs.recode_text))
    )

    ldl_meas = ldl_meas.annotate(**{
        drug: hl.cond(
            hl.is_defined(ldl_meas.diff_map),
            ldl_meas.diff_map.get(drug),
            hl.missing(hl.tint32)
        )
        for drug in unique_drugs
    })

    sources = ['GP', 'ASCEN']
    ldl_meas = ldl_meas.annotate(**{
        s: hl.int(ldl_meas.source.contains(s)) for s in sources
    })

    ldl_meas = ldl_meas.drop(ldl_meas.diff_map, ldl_meas.source)

    tb_name = "ldl-for-export.ht"
    url = f"dnax://{db_uri}/{tb_name}"
    ldl_meas.write(url, overwrite=True)
    ldl_meas = hl.read_table(url)

    return ldl_meas

