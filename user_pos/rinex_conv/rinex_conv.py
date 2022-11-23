import georinex as gr
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

def obs_read(obs_file, psuedo_file):
    obs = gr.load(obs_file)

    with open(psuedo_file, 'w') as f:
        #csv header
        f.write("GPS_Time,G02,G05,G10,G13,G15,G16,G18,G23,G24,G26,G27,G29,G32\n")

        for i in range(len(obs.time)):
            #make time into string
            time = np.datetime_as_string(obs.time[i])
            g2 = str(obs.sel(sv='G02')['L1'].values[i])
            g5 = str(obs.sel(sv='G05')['L1'].values[i])
            g10 = str(obs.sel(sv='G10')['L1'].values[i])
            g13 = str(obs.sel(sv='G13')['L1'].values[i])
            g15 = str(obs.sel(sv='G15')['L1'].values[i])
            g16 = str(obs.sel(sv='G16')['L1'].values[i])
            g18 = str(obs.sel(sv='G18')['L1'].values[i])
            g23 = str(obs.sel(sv='G23')['L1'].values[i])
            g26 = str(obs.sel(sv='G26')['L1'].values[i])
            g27 = str(obs.sel(sv='G27')['L1'].values[i])
            g29 = str(obs.sel(sv='G29')['L1'].values[i])
            g32= str(obs.sel(sv='G32')['L1'].values[i])

            #write data to file as csv
            f.write(time + ',' + g2 + ', ' + g5 + ',' + g10 + ',' + g13 + ',')
            f.write(g15 + ', ' + g16 + ',' + g18 + ',' + g23 + ',' + g26 + ',')
            f.write(g27 + ', ' + g29 + ',' + g32 + '\n')

def nav_read(nav_file):
    nav = gr.load(nav_file)

    for i in range(32):
        if i == 27:
            continue

        file_name = "mlfp325_csv/G" + str(i + 1) + "_nav.csv"
        sat = ""
        if(i < 9):
            sat = "G0" + str(i + 1)
        else:
            sat = "G" + str(i + 1)

        with open(file_name, 'w') as f:
            f.write('time,SVClockBias,SVClockDrift,SVclockDriftRate,IODE,Crs,DeltaN,M0,Cuc,Ecentricity,Cus,sqrtA,Toe,Cic,Io,omega,OmegaDot,IDOT,CodesL2,GPSWeek,L2Pflag,health,TGD,IODC,TransTime,FitIntvl\n')
            for i in range(len(nav.time)):
                time = np.datetime_as_string(nav.time[i])
                clockBias = str(nav.sel(sv=sat)['SVclockBias'].values[i])
                clockdrift = str(nav.sel(sv=sat)['SVclockDrift'].values[i])
                clockdriftrate = str(nav.sel(sv=sat)['SVclockDriftRate'].values[i])
                crs = str(nav.sel(sv=sat)['Crs'].values[i])
                deltan = str(nav.sel(sv=sat)['DeltaN'].values[i])
                m0 = str(nav.sel(sv=sat)['M0'].values[i])
                cuc = str(nav.sel(sv=sat)['Cuc'].values[i])
                ecentricity = str(nav.sel(sv=sat)['Eccentricity'].values[i])
                cus = str(nav.sel(sv=sat)['Cus'].values[i])
                sqrta = str(nav.sel(sv=sat)['sqrtA'].values[i])
                omega = str(nav.sel(sv=sat)['omega'].values[i])
                omegadot = str(nav.sel(sv=sat)['OmegaDot'].values[i])
                idot = str(nav.sel(sv=sat)['IDOT'].values[i])
                CodesL2 = str(nav.sel(sv=sat)['CodesL2'].values[i])
                GPSWeek = str(nav.sel(sv=sat)['GPSWeek'].values[i])
                L2Pflag = str(nav.sel(sv=sat)['L2Pflag'].values[i])
                health = str(nav.sel(sv=sat)['health'].values[i])
                TGD = str(nav.sel(sv=sat)['TGD'].values[i])
                IODC = str(nav.sel(sv=sat)['IODC'].values[i])
                TransTime = str(nav.sel(sv=sat)['TransTime'].values[i])
                FitIntvl = str(nav.sel(sv=sat)['FitIntvl'].values[i])

                f.write(time + ',' + clockBias + ', ' + clockdrift + ',' + clockdriftrate + ',' + crs + ',')
                f.write(deltan + ', ' + m0 + ',' + cuc + ',' + ecentricity + ',' + cus + ',')
                f.write(sqrta + ', ' + omega + ',' + omegadot + ',' + idot + ',' + CodesL2 + ',')
                f.write(GPSWeek + ', ' + L2Pflag + ',' + health + ',' + TGD + ',')
                f.write(IODC + ', ' + TransTime + ',' + FitIntvl + '\n')


def orbit_read(orbit_file):
    orbit = gr.load(orbit_file)

    for i in range(32):
        if i == 27:
            continue

        file_name = "mlfp325_csv/G" + str(i + 1) + "_orbit.csv"
        sat = ""
        if(i < 9):
            sat = "G0" + str(i + 1)
        else:
            sat = "G" + str(i + 1)

        with open(file_name, 'w') as f:
            f.write('time,X,Y,Z,clock,Vx,Vy,Vz,dclock\n')
            for i in range(len(orbit.time)):
                time = np.datetime_as_string(orbit.time[i])
                X = str(orbit.sel(sv=sat, ECEF="x")['position'].values[i])
                Y = str(orbit.sel(sv=sat, ECEF="y")['position'].values[i])
                Z = str(orbit.sel(sv=sat, ECEF="z")['position'].values[i])
                clock = str(orbit.sel(sv=sat)['clock'].values[i])
                Vx = str(orbit.sel(sv=sat, ECEF="x")['velocity'].values[i])
                Vy = str(orbit.sel(sv=sat, ECEF="y")['velocity'].values[i])
                Vz = str(orbit.sel(sv=sat, ECEF="z")['velocity'].values[i])
                dclock = str(orbit.sel(sv=sat)['dclock'].values[i])

                f.write(time + ',' + X + ', ' + Y + ',' + Z + ',' + clock + ',')
                f.write(Vx + ', ' + Vy + ',' + Vz + ',' + dclock + '\n')

def main():
    obs_read("mlfp325/mlfp3250.22o", "mlfp325_csv/obs.csv")
    nav_read("mlfp325/mlfp3250.22n")
    orbit_read("mlfp325/igu22371_180.sp3")

if __name__ == '__main__':
    main()
