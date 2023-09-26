from datetime import datetime


def merge_season(date, hemi, osi405yr=False):
    '''The winter and summer seasons are absolute, and the autumn
    and spring seasons are transition seasons. This returns the
    season and the season limits for any date.'''

    year = datetime.strftime(date, '%Y')
    if (hemi == 'nh' and
        date <= datetime.strptime('{}0501'.format(year), '%Y%m%d')):
        year = str(int(year) - 1)
    if (hemi == 'sh' and
        date <= datetime.strptime('{}0301'.format(year), '%Y%m%d')):
        year = str(int(year) - 1)

    mp = {}
    mp['nh'] = {}
    mp['nh']['winter'] = (datetime.strptime('{}1101'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0501'.format(str(int(year)+1)),
                                            '%Y%m%d'))
    mp['nh']['spring'] = (datetime.strptime('{}0501'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0601'.format(year), '%Y%m%d'))
    mp['nh']['summer'] = (datetime.strptime('{}0601'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1001'.format(year), '%Y%m%d'))
    mp['nh']['autumn'] = (datetime.strptime('{}1001'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1101'.format(year), '%Y%m%d'))
    mp['sh'] = {}
    mp['sh']['winter'] = (datetime.strptime('{}0401'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1001'.format(year), '%Y%m%d'))
    mp['sh']['spring'] = (datetime.strptime('{}1001'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1101'.format(year), '%Y%m%d'))
    mp['sh']['summer'] = (datetime.strptime('{}1101'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0301'.format(str(int(year)+1)),
                                            '%Y%m%d'))
    mp['sh']['autumn'] = (datetime.strptime('{}0301'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0401'.format(year), '%Y%m%d'))

    oy = {}
    oy['nh'] = {}
    oy['nh']['winter'] = (datetime.strptime('{}1001'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0401'.format(str(int(year)+1)),
                                            '%Y%m%d'))
    oy['nh']['spring'] = (datetime.strptime('{}0401'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0501'.format(year), '%Y%m%d'))
    oy['nh']['summer'] = (datetime.strptime('{}0501'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0901'.format(year), '%Y%m%d'))
    oy['nh']['autumn'] = (datetime.strptime('{}0901'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1001'.format(year), '%Y%m%d'))
    oy['sh'] = {}
    oy['sh']['winter'] = (datetime.strptime('{}0401'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1001'.format(year), '%Y%m%d'))
    oy['sh']['spring'] = (datetime.strptime('{}1001'.format(year), '%Y%m%d'),
                          datetime.strptime('{}1101'.format(year), '%Y%m%d'))
    oy['sh']['summer'] = (datetime.strptime('{}1101'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0301'.format(str(int(year)+1)),
                                            '%Y%m%d'))
    oy['sh']['autumn'] = (datetime.strptime('{}0301'.format(year), '%Y%m%d'),
                          datetime.strptime('{}0401'.format(year), '%Y%m%d'))

    if not osi405yr:
        merge_pers = mp[hemi]
    else:
        merge_pers = oy[hemi]
    for seas in ['winter', 'summer', 'spring', 'autumn']:
        if date >= merge_pers[seas][0] and date <= merge_pers[seas][1]:
            season = seas
            break

    return season, merge_pers[season]


def merge_frac(date, hemi, osi405yr=False):
    '''This takes one date and returns a fraction winter and a fraction
    summer for the merge'''

    season, lims = merge_season(date, hemi, osi405yr=osi405yr)
    if season == 'winter':
        return (1, 0)
    elif season == 'summer':
        return (0, 1)
    else:
        llim, ulim = lims
        lfrac = (date - llim) / (ulim - llim)
        ufrac = (ulim - date) / (ulim - llim)

        # Returns (winter fraction, summer fraction)
        if season == 'spring':
            return (ufrac, lfrac)
        elif season == 'autumn':
            return (lfrac, ufrac)
