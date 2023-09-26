import numpy as np
from datetime import datetime, date, time

def access_inst_field(I,f):
    content = I
    for ff in f.split('.'):
        try:
            content=content.__dict__[ff]
        except KeyError as k:
            raise KeyError("Instance object has no field named %s" % f)
        except Exception as e:
            raise Exception("Is not an Instance object (%s)" % e)
    return content

def access_dict_field(D,f):
    content = D
    for ff in f.split('.'):
        try:
            content=content[ff]
        except KeyError as k:
            raise KeyError("Dict object has no field named %s" % f)
        except Exception as e:
            raise Exception("Is not a Dict object (%e)" % e)
    return content

class selector:
    def __init__(self,neg=False):
        if not isinstance(neg,bool):
            raise TypeError("The neg= field must be bool (got {})".format(type(neg)))
        self.neg = neg

    def __str__(self):
        negstr = ''
        if self.neg:
            negstr='NOT'
        return negstr

    def select(self,M):
        if self.neg:
            return ~self.selectPos(M)
        else:
            return self.selectPos(M)

    def mask(self, M):
        if self.neg:
            return self.selectPos(M)
        else:
            return ~self.selectPos(M)

class single_selector(selector):
    def __init__(self, f, neg=False):
        selector.__init__(self,neg=neg)
        self.f   = f

    def __str__(self):
        return selector.__str__(self) + ' <{}>'.format(self.f)

    def getF(self, M):
        if len(self.f) == 0:
            ret = M
        else:
            try:
                ret = access_dict_field(M,self.f)
            except KeyError as k:
                raise KeyError("Dict object does not have field %s" % k) # forward KeyError
            except Exception as e:
                try:
                    ret = access_inst_field(M,self.f)
                except KeyError as k:
                    raise KeyError("Instance object does not have field %s" % k) # forward KeyError
                except Exception as e:
                    raise Exception("Attempt to apply a selector on field " +\
                                    "<%s> of an object that has no fields!" %\
                                    self.f);
        import scipy as sp
        if (type(ret) == sp.ma.core.MaskedArray):
            return ret.data
        else:
            return ret


# ###############################
# MULTI PURPOSE VALUE SELECTOR
# ###############################
class value_selector(single_selector):
    def __init__(self,v,f,neg=False,op='=='):
        single_selector.__init__(self,f,neg=neg)
        self.v = v
        self.op = op
    def __str__(self):
        return "%s VAL with val %s %s (%s)" % (single_selector.__str__(self), self.op,self.v, type(self.v))
    def selectPos(self, V):
        if self.op == '==':
            return self.getF(V) == self.v
        elif self.op == '>':
            return self.getF(V) > self.v
        elif self.op == '<':
            return self.getF(V) < self.v
        else:
            raise NotImplementedError('Operator {} is not implemented'.format(self.op))


# ###############################
# TEMPORAL SELECTORS
# ###############################
class time_selector(single_selector):
    def __init__(self, t, f, neg=False):
        single_selector.__init__(self,f,neg=neg)
        if isinstance(t,datetime):
            self.t = t
        elif isinstance(t,date):
            self.t = datetime.combine(t,time(0))
        else:
            raise TypeError("Try to initialize a time_selector with a {}: {}".format(type(t),t))
    def __str__(self):
        return "%s TS with timestamp %s" % (single_selector.__str__(self), self.t.isoformat())

class after_time_selector(time_selector):
    def __str__(self):
        return time_selector.__str__(self)+" AFTER"
    def selectPos(self, T):
        ret = self.getF(T) >= self.t
        #print self.getF(T)[~ret]
        return ret

class before_time_selector(time_selector):
    def __str__(self):
        return time_selector.__str__(self)+" BEFORE"
    def selectPos(self, T):
        ret = self.getF(T) <= self.t
        #print self.getF(T)[~ret]
        return ret

class period_selector(single_selector):
    def __init__(self, t0, t1, f, neg=False):
        single_selector.__init__(self,f,neg=neg)
        self.ts0 = after_time_selector(t0,f)
        self.ts1 = before_time_selector(t1,f)
    def __str__(self):
        return "%s PERIOD TS from %s to %s" % (single_selector.__str__(self), self.ts0.t, self.ts1.t)
    def __enter__(self):
        if self.ts1.t <= self.ts0.t:
            raise ValueError("Period is wrong: start time (%s) comes after stop time (%s)" % \
                            (self.ts0.t, self.ts1.t))
    def __exit__(self, exc_type, exc_value, exc_traceback):
        pass
    def selectPos(self, T):
        with self:
            r0 = (self.ts0).selectPos(T)
            r1 = (self.ts1).selectPos(T)
            r = r0 * r1
            return r

# ###############################
# TEXT/LABEL SELECTORS
# ###############################
class label_selector(single_selector):
    def __init__(self, s, f, neg=False):
        single_selector.__init__(self,f,neg=neg)
        self.s = s
    def __enter__(self):
        if self.s == '':
            raise ValueError("Label is zero-length string")
    def __exit__(self, exc_type, exc_value, exc_traceback):
        pass
    def __str__(self):
        return "%s SS with label <%s>" % (single_selector.__str__(self), self.s)

class exact_label_selector(label_selector):
    def __str__(self):
        return "%s EXACT MATCH" % label_selector.__str__(self)
    def selectPos(self, L):
        with self:
            #print(self.s,self.getF(L),type(self.s),type(self.getF(L)))
            return self.s == self.getF(L)

class nocase_label_selector(label_selector):
    def __str__(self):
        return "%s CASE-Insensitive MATCH" % label_selector.__str__(self)
    def selectPos(self, L):
        with self:
            up = self.s.upper()
            return np.array([up == l.upper() for l in self.getF(L)])

# ###############################
# GEOGRAPHICAL SELECTORS
# ###############################
class lat_selector(single_selector):
    def __init__(self, bbox, f, neg=False):
        single_selector.__init__(self,f,neg=neg)
        self.southernmost_lat = bbox[0]
        self.northernmost_lat = bbox[1]
    def __enter__(self):
        if self.southernmost_lat < -90. or self.southernmost_lat > 90.:
            raise ValueError("Latitudes must be in [-90:90] (got {})".format(self.southernmost_lat))
        if self.northernmost_lat < -90. or self.northernmost_lat > 90.:
            raise ValueError("Latitudes must be in [-90:90] (got {})".format(self.northernmost_lat))
        if self.southernmost_lat >= self.northernmost_lat:
            raise ValueError("The Southernmost lat must be lower than Northernmost one (got {} and {})".format(
                self.southernmost_lat,self.northernmost_lat))
    def __exit__(self, exc_type, exc_value, exc_traceback):
        pass
    def __str__(self):
        return "%s LAT with box: (%s < lat < %s)" % \
                (single_selector.__str__(self), self.southernmost_lat, self.northernmost_lat)
    def selectPos(self, P):
        with self:
            ret = (self.getF(P) >= self.southernmost_lat) * (self.getF(P) <= self.northernmost_lat)
            #print type(self.getF(P)), type(self.southernmost_lat), type(self.northernmost_lat)
            #print ret, type(ret)
            return (self.getF(P) >= self.southernmost_lat) * (self.getF(P) <= self.northernmost_lat)

class lon_selector(single_selector):
    def __init__(self, bbox, f, neg=False):
        single_selector.__init__(self,f,neg=neg)
        # force all lons in [-180:+180[
        self.westernmost_lon  = (bbox[0] + 180) % (360) - 180
        self.easternmost_lon  = (bbox[1] + 180) % (360) - 180
    def __enter__(self):
        if self.westernmost_lon >= self.easternmost_lon:
            raise ValueError("The Westernmost lon must be lower than Eastermost one (got {} and {})".format(
                self.westernmost_lon,self.easternmost_lon))
    def __exit__(self, exc_type, exc_value, exc_traceback):
        pass
    def __str__(self):
        return "%s LON with box: (%s < lon < %s)" % \
                (single_selector.__str__(self), \
                 self.westernmost_lon, self.easternmost_lon)
    def selectPos(self, P):
        with self:
            return (self.getF(P) >= self.westernmost_lon) * (self.getF(P) <= self.easternmost_lon)

# ###############################
# COMBINING SELECTORS (AND, OR)
# ###############################
class combined_selector(selector):
    def __init__(self,s1,s2,op,neg=False):
        #if not isinstance(s1,single_selector) or not isinstance(s2,single_selector):
        #    raise ValueError("A combined selector can only be build from 2 matchup selectors")
        if not op.lower() in ('and','or',):
            raise ValueError("A combined selector can only have 'and' or 'or' operator")
        selector.__init__(self,neg=neg)
        self.s1 = s1
        self.s2 = s2
        self.op = op.lower()

    def __str__(self):
        return selector.__str__(self) + "[ ({}) {} ({}) ]".format(self.s1.__str__(),self.op.upper(),self.s2.__str__())

    def selectPos(self,W):
        if self.op == 'and':
            return self.s1.selectPos(W) * self.s2.selectPos(W)
        else:
            return self.s1.selectPos(W) + self.s2.selectPos(W)

if __name__ == '__main__':

    import scipy as sp
    from scipy import ma
    import numpy as np
    import os
    import base64
    from datetime import timedelta

    print("TIME SELECTOR")
    print('='*20)
    p_beg = datetime(2009,2,13)
    p_end = datetime(2009,2,17)

    ps = period_selector(p_beg,p_end,'t0')
    print("ps: %s" % ps)

    mT = {'t0': sp.array([p_beg - timedelta(hours=12) + timedelta(hours=h) for h in range(0,4*24+1)])}
    print("MT :", mT)

    print("Select: %s" % ps.select(mT))
    print("Mask: %s" % ps.mask(mT))

    print("LABEL SELECTOR")
    print('='*20)
    ls = exact_label_selector('itp21','id')
    print("ls: %s" % ls)

    mL = {'id': sp.array(['qbc', '23mi', 'itp21', 'ij56'])}
    print("ML: ", mL['id'])

    print("Select: %s" % ls.select(mL))
    print("Mask: %s" % ls.mask(mL))

    print("LABEL SELECTOR (NOT)")
    print('='*20)
    notls = exact_label_selector('23mi','id',neg=True)
    print("not-ls: %s" % notls)

    print("ML: ", mL['id'])

    print("Select: %s" % notls.select(mL))
    print("Mask: %s" % notls.mask(mL))

    print("VALUE SELECTOR")
    print('='*20)
    ms = value_selector(17,'flag')
    print(ms)

    mL = {'flag': sp.array([0,0,13,17,0,0,12,0,0,17])}
    print("Select: %s" % ms.select(mL))
    print("Mask: %s" % ms.mask(mL))


    print("GEO SELECTOR")
    print('='*20)
    latS = lat_selector((70.,85.),'lat')
    print(latS)
    lonS = lon_selector((-30.,10.),'lon')
    print(lonS)

    pos = {'lon': sp.arange(-45,45,10.)}
    pos['lat'] = sp.array([81.]*len(pos['lon']))
    print(pos)

    geoS = combined_selector(latS,lonS,'and',neg=True)
    print(geoS)
    print("Select: %s" % geoS.select(pos))
    print("Mask: %s" % geoS.mask(pos))
