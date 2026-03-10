import pandas as pd

input_excel_file = '../data/lkps/all_lkps_maps_v4.xlsx'

all_sheets = pd.read_excel(input_excel_file, sheet_name=None)
sheets = all_sheets.keys()
print(sheets)

for sheet in sheets:
    csv_file_name = f"{sheet.replace(' ', '_').replace('/', '_')}.csv"
    data = all_sheets[sheet]
    data.to_csv(csv_file_name, index=False)
    print(f"Zapisano arkusz '{sheet}' do pliku: {csv_file_name}")
