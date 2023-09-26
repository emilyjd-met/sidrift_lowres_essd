'''

NAME
       find_matchups_files.py

SYNOPSIS
       python find_matchups_files [OPTIONS]

DESCRIPTION
       .......


OPTIONS
       A summary of the options supported by find_matchups_files.py
       is included below.

       -dato      Selects the period given between first_date
                  and last_date.
                  ex: python find_matchups_files -dato,2012,1,2013,12

       -H1        Selects the 6 first months of the given year
       -H2        Selects the 6 last months of the given year
                  ex: python find_matchups_files -H2,2013

       -Q1        Selects the first quarter of the given year
       -Q2        Selects the second quarter of the given year
       -Q3        Selecta the third quarter of the given year
       -Q4        Selects the last quarter of the given year
                  ex: python find_matchups_files -Q3,2013

       -Y         Selects the given year
                  ex: python find_matchups_files -Y,2013

       -M1        Selects the first month of the given year
        .
        .
        .
       -M12       Selects the last month of the given year
                  ex: python find_matchups_files -M12,2013

'''

import sys
import os
import glob
import datetime
from dateutil import rrule, relativedelta

def merge_periods(times):
    """
       Re-order and merge successive time periods

       Test merging of a single time period
       >>> s1 = datetime.date(2013,1,1)
       >>> e1 = datetime.date(2013,1,31)
       >>> arg = [[s1,e1],]
       >>> merge_periods(arg) == arg
       True
       >>> arg_tuple = [(s1,e1),]
       >>> merge_periods(arg_tuple)
       [[datetime.date(2013, 1, 1), datetime.date(2013, 1, 31)]]
       >>> arg_tuple = ((s1,e1),)
       >>> merge_periods(arg_tuple)
       [[datetime.date(2013, 1, 1), datetime.date(2013, 1, 31)]]

       Test merging of two consecutive time periods into one
       >>> s1 = datetime.date(2013,1,1)
       >>> e1 = datetime.date(2013,1,31)
       >>> s2 = datetime.date(2013,2,1)
       >>> e2 = datetime.date(2013,2,28)
       >>> arg = [[s1,e1],[s2,e2]]
       >>> merge_periods(arg)
       [[datetime.date(2013, 1, 1), datetime.date(2013, 2, 28)]]

       Test merging of three time periods into two
       >>> s1 = datetime.date(2013,1,1)
       >>> e1 = datetime.date(2013,1,31)
       >>> s2 = datetime.date(2013,2,1)
       >>> e2 = datetime.date(2013,2,28)
       >>> s3 = datetime.date(2013,4,5)
       >>> e3 = datetime.date(2013,4,25)
       >>> arg = [[s1,e1],[s2,e2],[s3,e3]]
       >>> merge_periods(arg)
       [[datetime.date(2013, 1, 1), datetime.date(2013, 2, 28)], [datetime.date(2013, 4, 5), datetime.date(2013, 4, 25)]]

       Test re-ordering of periods
       >>> s1 = datetime.date(2013,1,1)
       >>> e1 = datetime.date(2013,1,31)
       >>> s2 = datetime.date(2012,7,1)
       >>> e2 = datetime.date(2012,9,15)
       >>> arg1 = [[s1,e1],[s2,e2]]
       >>> arg2 = [[s2,e2],[s1,e1]]
       >>> merge_periods(arg1) == merge_periods(arg2)
       True

       Test re-ordering AND merging of time periods
       >>> s1 = datetime.date(2013,1,1)
       >>> e1 = datetime.date(2013,1,31)
       >>> s2 = datetime.date(2012,7,1)
       >>> e2 = datetime.date(2012,12,31)
       >>> arg1 = [[s1,e1],[s2,e2]]
       >>> arg2 = [[s2,e2],[s1,e1]]
       >>> merge_periods(arg1) == merge_periods(arg2)
       True
       >>> len(merge_periods(arg1)) == 1
       True

    """
    periods = []
    times = list(times)
    times = sorted(times)
    saved = list(times[0])
    for st, en in sorted([sorted(t) for t in times]):
        if st <= (saved[1]+datetime.timedelta(days=1)):
            saved[1] = max(saved[1], en)
        else:
            periods.append(saved[:])
            saved[0] = st
            saved[1] = en

    periods.append(saved[:])

    return periods

def decode_period(period):
    """
       decode a set of periods separated by + signs. Also
       sorts and merge together consecutive (and overlapping) periods

       Test a single period
       >>> p1 = 'H2,2013'
       >>> decode_period(p1)
       ([datetime.date(2013, 7, 1)], [datetime.date(2013, 12, 31)])

       Test two distinct single periods
       >>> p1 = 'H2,2013+M2,2014'
       >>> d1,d2 = decode_period(p1)
       >>> d1
       [datetime.date(2013, 7, 1), datetime.date(2014, 2, 1)]
       >>> d2
       [datetime.date(2013, 12, 31), datetime.date(2014, 2, 28)]

       Test two successive single periods
       >>> p1 = 'H2,2013+M1,2014+H2,2014'
       >>> d1,d2 = decode_period(p1)
       >>> d1
       [datetime.date(2013, 7, 1), datetime.date(2014, 7, 1)]
       >>> d2
       [datetime.date(2014, 1, 31), datetime.date(2014, 12, 31)]
    """

    # decode strings
    date1 = []
    date2 = []
    for p in period.split('+'):
        d1, d2 = decode_one_period(p)
        date1.append(d1)
        date2.append(d2)

    # re-order and merge successive time periods
    args = zip(date1,date2)
    mper = merge_periods(zip(date1,date2))
    merged_date1 = []
    merged_date2 = []
    for i in range(len(mper)):
        merged_date1.append(mper[i][0])
        merged_date2.append(mper[i][1])

    return merged_date1, merged_date2

def decode_one_period(period):
    """
       return date_start and date_stop from a period string

       Test Half-Year
       >>> decode_one_period('H2,2013')
       (datetime.date(2013, 7, 1), datetime.date(2013, 12, 31))
       >>> decode_one_period('H1,2012')
       (datetime.date(2012, 1, 1), datetime.date(2012, 6, 30))

       Test Quarters
       >>> decode_one_period('Q1,2012')
       (datetime.date(2012, 1, 1), datetime.date(2012, 3, 31))
       >>> decode_one_period('Q2,2012')
       (datetime.date(2012, 4, 1), datetime.date(2012, 6, 30))
       >>> decode_one_period('Q3,2012')
       (datetime.date(2012, 7, 1), datetime.date(2012, 9, 30))
       >>> decode_one_period('Q4,2012')
       (datetime.date(2012, 10, 1), datetime.date(2012, 12, 31))

       Test Year
       >>> decode_one_period('Y,2009')
       (datetime.date(2009, 1, 1), datetime.date(2009, 12, 31))

       Test Month
       >>> decode_one_period('M2,2008')
       (datetime.date(2008, 2, 1), datetime.date(2008, 2, 29))
       >>> decode_one_period('M2,2007')
       (datetime.date(2007, 2, 1), datetime.date(2007, 2, 28))
       >>> decode_one_period('M7,2010')
       (datetime.date(2010, 7, 1), datetime.date(2010, 7, 31))
       >>> decode_one_period('M11,2019')
       (datetime.date(2019, 11, 1), datetime.date(2019, 11, 30))

   """

    mm = []
    yy = []
    mm1 = []
    mm2 = []

    if len(period) == 0:
        raise ValueError("Cannot decode empty period string")

    i1 = period.split(',')
    i2 = int(period.split(',')[1])
    if i1[0] == 'dato':
        if len(i1) != 5:
            raise ValueError("wrong OPTION")
        yy = range(int(period.split(',')[1]),int(period.split(',')[3])+1)
        mm1 = int(period.split(',')[2])
        mm2 = int(period.split(',')[4])
    elif i1[0][0] == 'Q':
        if len(i1) != 2:
            raise ValueError("wrong OPTION")
        Q = int(i1[0][1:])
        yy = i2
        for m in range(1, 13):
            if (m-1)//3+1 == Q:
                mm.append(m)
        mm1 = min(mm)
        mm2 = max(mm)
    elif i1[0][0] == 'H':
        if len(i1) != 2:
            raise ValueError("wrong OPTION")
        H = i1[0][1:]
        yy = i2
        if H == '1':
            mm1 = 1
            mm2 = 6
        elif H == '2':
            mm1 = 7
            mm2 = 12
    elif i1[0][0] == 'Y':
        if len(i1) != 2:
            raise ValueError("wrong OPTION")
        yy = i2
        mm1 = 1
        mm2 = 12
    elif i1[0][0] == 'M':
        if len(i1) != 2:
            raise ValueError("wrong OPTION")
        yy = i2
        mm1 = int(i1[0][1:])
        mm2 = mm1
    else:
        raise ValueError("This OPTIONS is not possible")

    import calendar
    if type(yy) == list:
        dato1=datetime.date(min(yy),mm1,1)
        dato2=datetime.date(max(yy),mm2,calendar.monthrange(max(yy),mm2)[1])
    elif type(yy) == int:
        dato1=datetime.date(yy,mm1,1)
        dato2=datetime.date(yy,mm2,calendar.monthrange(yy,mm2)[1])
    return(dato1,dato2)

def get_previous_period(period):
    """
    >>> get_previous_period('Y,2010')
    'Y,2009'
    >>> get_previous_period('M11,2019')
    'M10,2019'
    >>> get_previous_period('M1,2008')
    'M12,2007'
    >>> get_previous_period('H1,2005')
    'H2,2004'
    >>> get_previous_period('H2,2005')
    'H1,2005'
    >>> get_previous_period('Q1,2005')
    'Q4,2004'
    >>> get_previous_period('Q3,2010')
    'Q2,2010'
    """
    pStr,pYear = period.split(',')
    pCode  = pStr[0]

    # special case for Yearly period
    if pCode == 'Y':
        return 'Y,{}'.format(int(pYear)-1)

    # general case for M, Q, H
    pValue = int(pStr[1:])
    if pValue < 1:
        raise ValueError("Invalid period string {}".format(period,))
    if pCode == 'H':
        maxValue = 2
    elif pCode == 'Q':
        maxValue = 4
    elif pCode == 'M':
        maxValue = 12
    else:
        raise ValueError("Invalid period string {}".format(period,))

    newValue = pValue - 1
    newYear  = int(pYear)
    if newValue < 1:
        newValue = maxValue
        newYear -= 1

    return '{}{},{}'.format(pCode,newValue,newYear)

def split_period_in_months(date1,date2,reverse=False):
    """
        Split a time period (date1,date2) in monthly period:

        Note: the reverse keyword does allow to loop through the months
           backward in time (but the months are still defined from first to last day)

        >>> split_period_in_months([datetime.date(2012,2,15)],[datetime.date(2012,2,16)])
        [[datetime.date(2012, 2, 1), datetime.date(2012, 2, 29)]]
        >>> split_period_in_months([datetime.date(2010,2,15)],[datetime.date(2010,3,3)])
        [[datetime.date(2010, 2, 1), datetime.date(2010, 2, 28)], [datetime.date(2010, 3, 1), datetime.date(2010, 3, 31)]]
        >>> split_period_in_months([datetime.date(2010,2,15)],[datetime.date(2010,3,3)],reverse=True)
        [[datetime.date(2010, 3, 1), datetime.date(2010, 3, 31)], [datetime.date(2010, 2, 1), datetime.date(2010, 2, 28)]]

    """

    # Select the start and end dates from a list if a list is given
    if isinstance(date1, list):
        date1 = date1[0]
    if isinstance(date2, list):
        date2 = date2[-1]

    sub_periods = []
    for p in range(len(date1)):
        fday_of_fmon = datetime.date(date1[p].year,date1[p].month,1)
        for dt in rrule.rrule(rrule.MONTHLY, dtstart=fday_of_fmon, until=date2[p]):
            sub_periods.append([dt.date(),(dt+relativedelta.relativedelta(months=1)-datetime.timedelta(days=1)).date()])

    if reverse:
        sub_periods = sub_periods[::-1]

    return sub_periods


def split_period_in_years(date1, date2, startdate=(1, 1)):

    # The startdate is a tuple of (month, year) that the year runs from

    # Select the start and end dates from a list if a list is given
    if isinstance(date1, list):
        date1 = date1[0]
    if isinstance(date2, list):
        date2 = date2[-1]

    sub_periods = []
    for p in range(len(date1)):
        fyear = datetime.date(date1[p].year, startdate[0], startdate[1])
        for dt in rrule.rrule(rrule.YEARLY, dtstart=fyear, until=date2[p]):
            sub_periods.append([dt.date(),
                                (dt + relativedelta.relativedelta(years=1)
                                 - datetime.timedelta(days=1)).date()])

    return sub_periods


def iterate_daily_in_period(Lstart,Lstop):
    """
       Yield daily increments in a period.
    """

    if len(Lstart) != len(Lstop):
        raise ValueError("The Lstart and Lstop arrays must have the same length")

    for p in range(len(Lstart)):
        p_b = Lstart[p]
        p_e = Lstop[p]
        delta = p_e - p_b
        for i in range(delta.days + 1):
            loop_date = p_b + datetime.timedelta(days=i)
            yield loop_date


def find_matchup_files(Lstart,Lstop,datadir='/disk1/osisaf/mdb-icedrift',instr='multi-oi',
                       lev='lev3',coll='NN3D',output=None,area='nh',
                       reproc_instr=False,reproc_date=False,
                       channels='',swaths=False,tspan=2, yearly=False,
                       verbose=False):
    """ returns list of matchup files """

    instr_ch = instr
    if instr == 'wind':
        instr_ch = 'ecmwf-era5_none_wind'
    elif instr != 'multi-oi':
        instr_ch += '_' + channels
    if instr == 'wind':
        simplex = ''
    else:
        simplex = '_simplex'

    timepatt = '12'
    if swaths:
        timepatt = '??????'

    list_files = []
    for ip in iterate_daily_in_period(Lstart,Lstop):
        ip2 = ip + datetime.timedelta(days = tspan)

        year  = ip2.year
        month = ip2.month
        day   = ip2.day

        if yearly:
            name = 'matchup_{}_{}_{:%Y}_{}col.nc'.format(
                instr, area, ip, coll)
        else:
            name = 'matchup_{}{}_{}_{}_{:%Y%m%d}{}w_{:%Y%m%d}{}w_{}col.nc'.format(instr_ch, simplex, lev, area, ip, timepatt, ip2, timepatt, coll)

        # subdirs to look for the matchup files
        subdirs = ['.',]

        if reproc_instr:
            subdirs = [instr,]

        if reproc_date:
            subdirs = [subdirs[0]+'/{d:%Y}/{d:%m}/'.format(d=ip2),\
                       subdirs[0]+'/{d:%Y}/{d:%m}/{d:%d}/'.format(d=ip2),]

        path = ["/".join([datadir,sd,name]) for sd in subdirs]
        for p in path:
            mf = glob.glob( p )
            if verbose and len(mf) > 0:
                print("Found {} files with {}".format(len(mf),p))


            list_files.extend(mf)

    # unique the list of files
    u_list_files = list(set(list_files))

    # dump to file (for debugging)
    if output is not None:
        fout = open(output, 'w')
        for string in u_list_files:
            fout.write('{0:s}\n'.format(string))
        fout.close()
    else:
        #print("Do not write list of files to file")
        pass

    if verbose:
        print("Found %d matchup files" % (len(u_list_files),))

    return u_list_files

def main():
    """
        Main loop.
    """
    if len(sys.argv) == 1:
        print("ERROR: You must provide a period string as parameter")
        sys.exit(1)
    dato1,dato2 = decode_period(str(sys.argv[1]))
    print(dato1,dato2)
    l = find_matchup_files(dato1,dato2)
    print(l)

if (__name__ == '__main__'):
    import doctest
    doctest.testmod()

    main()
