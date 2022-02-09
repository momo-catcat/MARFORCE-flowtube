# basic input for tubes
# MION inlet
# 09.10 the frist SA cali with 41 cm for 3/4 inch and 58.5 cm for 1 inch first tower Br API9
# 10.28 the second SA cali with 50 cm for 3/4 inch and 68 cm for 1 inch we have Y pieces first tower Br API9
# 11.18 the third SA cali with 50 cm for 3/4 inch and 66 cm for 1 inch we have Y pieces and second tower Br API9
# 01.04 the fourth SA cali with 10 cm for 3/4 inch and 78 cm for 1 inch we have Y pieces and first tower NO3 tower API9
# NO3 inlet
# 01.27 the fifth SA cali with 26cm for 3/4 inch NO3 inlet API9
# MION inlet
# 02.07 the sixth SA cali with 10 cm for 3/4 inch and 61 cm for 1 inch first tower Br karsa

def inputs_va(date):
    if date == '09.10':
        flag_tube = '2'
        R1 = 0.78
        L1 = 41
        R2 = 1.04
        L2 = 58.5
        file = 'H2O_1.csv'
        s1 = '1.csv'
        s2 = '1.txt'
    elif date == '10.28':
        flag_tube = '3'
        R1 = 0.78
        L1 = 50
        R2 = 1.04
        L2 = 68
        file = 'H2O_2.csv'
        s1 = '2.csv'
        s2 = '2.txt'
    elif date == '11.18':
        flag_tube = '3'
        R1 = 0.78
        L1 = 50
        R2 = 1.04
        L2 = 66
        file = 'H2O_3.csv'
        s1 = '3.csv'
        s2 = '3.txt'
    elif date == '01.04':
        flag_tube = '3'
        R1 = 0.78
        L1 = 10
        R2 = 1.04
        L2 = 78
        file = 'H2O_4.csv'
        s1 = '4.csv'
        s2 = '4.txt'
    elif date == '01.27':
        flag_tube = '1'
        R1 = 0.78
        L1 = 10
        R2 = 0
        L2 = 0
        file = 'H2O_1.csv'
        s1 = '5.csv'
        s2 = '5.txt'
    else:
        flag_tube = '2'
        R1 = 0.78
        L1 = 10
        R2 = 1.04
        L2 = 61
        file = 'H2O_1.csv'
        s1 = '6.csv'
        s2 = '6.txt'
    return R1,L1,R2,L2,flag_tube,file,s1,s2

