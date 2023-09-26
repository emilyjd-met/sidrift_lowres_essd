from datetime import datetime, timedelta
import calendar


def period_def(period, year):

    periods = {'summer': [datetime(year, 6, 1), datetime(year, 9, 30)],
               'winter': [datetime(year-1, 10, 1), datetime(year, 5, 31)],
               'jan': [datetime(year, 1, 1), datetime(year, 1, 31)],
               'feb': [datetime(year, 2, 1),
                       datetime(year, 2, 29) if calendar.isleap(year) else datetime(year, 2, 28)],
               'mar': [datetime(year, 3, 1), datetime(year, 3, 31)],
               'apr': [datetime(year, 4, 1), datetime(year, 4, 30)],
               'may': [datetime(year, 5, 1), datetime(year, 5, 31)],
               'jun': [datetime(year, 6, 1), datetime(year, 6, 30)],
               'jul': [datetime(year, 7, 1), datetime(year, 7, 31)],
               'aug': [datetime(year, 8, 1), datetime(year, 8, 31)],
               'sep': [datetime(year, 9, 1), datetime(year, 9, 30)],
               'oct': [datetime(year, 10, 1), datetime(year, 10, 31)],
               'nov': [datetime(year, 11, 1), datetime(year, 11, 30)],
               'dec': [datetime(year, 12, 1), datetime(year, 12, 31)]
       }

    if period not in periods.keys():
        raise(ValueError, "Unrecognised time period: {}\n".format(period))

    start_date = periods[period][0]
    end_date = periods[period][1]

    return start_date, end_date


def param_find(pdate):
    '''Find the two bordering parameter files with weightings'''

    # The midpoint of the given day is used (it is assumed that pdate
    # has no time field)
    pdate0 = pdate
    pdate = pdate.replace(hour=12)

    months = ['jan', 'feb', 'mar', 'apr', 'may', 'jun', 'jul', 'aug',
              'sep', 'oct', 'nov', 'dec']

    # Find which month the date is in and the offset to the mid-date of
    # the month
    for month in months:
        start_date, end_date = period_def(month, int(pdate.strftime('%Y')))
        if ((pdate0 >= start_date) and (pdate0 <= end_date)):
            mon1 = month
            break
    days1 = end_date - start_date + timedelta(days=1)
    mid_date1 = start_date + (0.5 * days1)
    offset1 = mid_date1 - pdate

    # Check whether it is the month before or after that is required,
    mon2_yr = int(pdate.strftime('%Y'))
    mon1_index = months.index(mon1)
    if offset1 >= timedelta(0):
        mon2_index = mon1_index - 1
        if mon2_index == -1:
            mon2_index = 11
    else:
        offset1 = offset1 * -1
        mon2_index = mon1_index + 1
        if mon2_index == 12:
            mon2_index = 0
            # Rolling over the year if necessary
            mon2_yr += 1
    mon2 = months[mon2_index]

    # Find the offset to the centre date of month 2
    start_date, end_date = period_def(mon2, mon2_yr)
    days2 = end_date - start_date + timedelta(days=1)
    mid_date2 = start_date + (0.5 * days2)
    offset2 = mid_date2 - pdate
    if offset2 < timedelta(0):
        offset2 = offset2 * -1

    weight1 = offset2.total_seconds() / timedelta(days=1).total_seconds()
    weight2 = offset1.total_seconds() / timedelta(days=1).total_seconds()

    return mon1, weight1, mon2, weight2


def test():

    p1 = datetime.strptime('20190214', '%Y%m%d')
    m1, w1, m2, w2 = param_find(p1)
    print("Date: {} * mon1: {} * weight1: {} * mon2: {} * weight2: {}".format(
        p1, m1, w1, m2, w2))

    p1 = datetime.strptime('20200215', '%Y%m%d')
    m1, w1, m2, w2 = param_find(p1)
    print("Date: {} * mon1: {} * weight1: {} * mon2: {} * weight2: {}".format(
        p1, m1, w1, m2, w2))

    p1 = datetime.strptime('20170305', '%Y%m%d')
    m1, w1, m2, w2 = param_find(p1)
    print("Date: {} * mon1: {} * weight1: {} * mon2: {} * weight2: {}".format(
        p1, m1, w1, m2, w2))

    p1 = datetime.strptime('20181219', '%Y%m%d')
    m1, w1, m2, w2 = param_find(p1)
    print("Date: {} * mon1: {} * weight1: {} * mon2: {} * weight2: {}".format(
        p1, m1, w1, m2, w2))

    p1 = datetime.strptime('20200530', '%Y%m%d')
    m1, w1, m2, w2 = param_find(p1)
    print("Date: {} * mon1: {} * weight1: {} * mon2: {} * weight2: {}".format(
        p1, m1, w1, m2, w2))


if __name__ == '__main__':

    test()
