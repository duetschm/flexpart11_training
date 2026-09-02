import os
import sys

from correlation_vis import correlation_plt


INPUT_FILES = ['./correlation.txt']
TUPLE_DOMAIN = None 
COL_SCALE = 'vast'
PLT_PROJ = 'cyl'
RUN_TYPE = 'single'
PLT_TYPE = 'conc_depo'
PLT_FIELD = 'corr'
PLT_MARKER = True
PLOT_TITLE = None
PLOT_TITLE_PART = None
NUM_TIME_STEP = 120
STATION_FILE = './station_file'
CORR_FILES = [input_file + '.srm' for input_file in INPUT_FILES]

# access inforamtion on additional markers
f = open('./markers.dat','r')
locations = f.readlines()
f.close()

# split file in srs and station
for input_file, srs_file in zip(INPUT_FILES, CORR_FILES):
    with open(input_file, 'r') as brioud_file:
        lines = brioud_file.readlines()
    # access information on involved IMS stations
    with open(STATION_FILE, 'w') as station_file:
        for line_idx, line in enumerate(lines):
            line_list = line.split()
            if len(line_list) < 3:
                if line_idx > 1:
                    line_list.append(' station' + str(line_idx - 1) + '\n')
                    station_file.write(' '.join(line_list))
                    if line_idx == 2:
                        lon_stat1 = float(line_list[0])
                        lat_stat1 = float(line_list[1])
                    del line_list
                else:
                    number = int(line)
                    number = number + len(locations)  
                    station_file.write(str(number) + '\n')
        for loc in locations:
            stat_lon_lat = loc.split()
            lon = stat_lon_lat[1]
            lat = stat_lon_lat[2]
            name = stat_lon_lat[0]
            station_file.write(lon + ' ' + lat + ' ' + name + '\n')
    lines_to_skip = int(lines[1])
    with open(srs_file, 'w') as srs_file_write:
        srs_file_write.write(lines[0])
        srs_file_write.writelines(lines[2 + lines_to_skip:])
STATION = {'lat': lat_stat1, 'lon': lon_stat1, 'station_id': ''}
stations = [STATION] * len(INPUT_FILES)
correlation_plt(
input_files=CORR_FILES,
stations=stations,
tuple_domain=TUPLE_DOMAIN,
col_scale=COL_SCALE,
plt_proj=PLT_PROJ,
run_type=RUN_TYPE,
plt_type=PLT_TYPE,
mark_stations=PLT_MARKER,
plot_title_part=PLOT_TITLE,
station_file=STATION_FILE,
plot_title=PLOT_TITLE_PART,
num_time_steps=NUM_TIME_STEP)

print('Plot psr is finished!')

