import sys
import os
import re
import scipy as sp
from scipy import ma
import numpy as np
import traceback
from operator import attrgetter
from netCDF4 import Dataset, num2date, chartostring

import matchups_selector as msel

def read_list_of_matchup_files(list_file, exit=True, sensor=None):
    """
       Read names of matchup files from file 'list_file', containing
       the full path to a macthup file per line.

       RETURN:
           the list of filenames

    """

    # open the file
    try:
        fh = open(list_file, "r")
    except Exception as e:
        raise Exception(e)

    # read it line-by-line
    lines = fh.readlines()

    # close
    fh.close()

    # check that each line has correct format
    all_ok = True
    nl = 1
    is_matchup_file_pattern = '^matchup_'
    if (sensor):
        is_matchup_file_pattern = is_matchup_file_pattern + sensor + '_'
    is_matchup_file_pattern = is_matchup_file_pattern + '.*' + '\.nc'
    is_matchup_file_re = re.compile(is_matchup_file_pattern)
    files = []
    for l in lines:
        ll = l.split()
        if (len(ll) != 1):
            raise ValueError("File %s not in valid format at line %s" \
                             % (list_file, nl))
        fname = ll[0]
        if not re.search(is_matchup_file_re, os.path.basename(fname)):
            msg = "Not a valid filename at line %s" % nl
            all_ok = False
            if exit:
                raise ValueError(msg)
            else:
                print("ERROR: %s" % msg)

        files.append(fname)
        nl = nl+1

    if not all_ok:
        raise ValueError("At least one filename is not valid in %s" %\
                         list_file)

    # return the list of files
    return files

class sid_collocation:
    def __init__(self):
        self.dist     = ma.array([])
        self.diff_t0  = ma.array([])
        self.diff_t1  = ma.array([])
        self.diff_dt  = ma.array([])
    def __str__(self):
        return "(VAL-PRD) dist:%s km diff_t0:%s diff_t1:%s diff_dt:%s" % \
                (self.dist, self.diff_t0, self.diff_t1, self.diff_dt)
    def _applymask(self, mask):

        for f in ('dist','diff_t0','diff_t1','diff_dt'):
            setattr(getattr(self, f),'mask',np.ma.nomask)
            getattr(self, f)[mask] = np.ma.masked


    def _compute(self, val, prd):
        # Compute collocation information
        from pyproj import Geod
        from datetime import timedelta
        g = Geod(ellps='sphere')
        self.dist     = ma.array(g.inv(val.lon0,val.lat0,prd.lon1,prd.lat1)[2])/1000
        self.diff_t0  = ma.array(val.t0 - prd.t0)
        self.diff_t1  = ma.array(val.t1 - prd.t1)
        self.diff_dt  = ma.array(val.dt - prd.dt)

class sid_estimates:

    fields = ('id','net','src','t0','t1','dt','lat0','lon0','dX','dY','sX','sY','cXY','lat1','lon1','dir','len','flg','relangle')

    def __init__(self):
        for f in self.fields:
            setattr(self, f, ma.array([]))

    def __str__(self):
        return "id:%s net:%s src:%s t0:%s t1:%s dt:%s lat0:%s lon0:%s dX:%s dY:%s lat1:%s lon1:%s len:%s dir:%s flg:%s relangle:%s" % \
                (self.id, self.net, self.src, self.t0, self.t1, self.dt,self.lat0, self.lon0, self.dX, self.dY,\
                self.lat1, self.lon1, self.len, self.dir, self.flg, self.relangle)

    def __enter__(self):
        """ Basic tests and validation of the object """
        self._compute_extra_sid_estimates_fields()
        # test that all arrays are mono-dimensional
        r_id   = self.id.shape
        for f in self.fields:
            if getattr(self,f).shape != r_id:
                raise ValueError("Not all matchup fields have same shape! e.g. {} {}".format(f,getattr(self,f).shape))

    def __exit__(self, exc_type, exc_value, exc_traceback):
        """ Shutdown nicely """
        pass

    def _applymask(self, mask):

        if mask is not ma.nomask:
            if not isinstance(mask,np.ndarray) or mask.dtype != bool:
                raise ValueError("_applymask() expects an ndarray of bools, got {} of {}".format(type(mask),mask.dtype))

        for f in self.fields:
            setattr(getattr(self, f),'mask',np.ma.nomask)
            getattr(self, f)[mask] = np.ma.masked


    def join(self,x):

        if isinstance(x,sid_estimates):
            x = (x,)

        for f in self.fields:
            all_f_fields = list(map(attrgetter(f),x))
            setattr(self, f, ma.concatenate((getattr(self, f), *all_f_fields)))

    def write_to_matchup_file(self, rootgrp, type):
        if (type is not 'val' and type is not 'prd'):
            raise ValueError("ERROR: 'type' (3rd param) must be 'val' or 'prd' in call to routine" % __name__)

        if 'matchup' not in rootgrp.dimensions():
            raise ValueError("ERROR: refuse to write <{}> because netCDF file is missing 'matchup'".format(type))

        lmatchup = len(rootgrp.dimensions()['matchup'])
        if self.lat0.count() != lmatchup:
            raise ValueError("ERROR: error, lengths do not match ({} != {})".format(self.lat0.count(), lmatchup))

    def load_from_matchup_file(self, rootgrp, type):
        if (type is not 'val' and type is not 'prd'):
            raise ValueError("ERROR: 'type' (3rd param) must be 'val' or 'prd' in call to routine" % __name__)

        try:
            self.lon0 = rootgrp.variables[type+'B_lon'][:].astype(float)
            self.lat0 = rootgrp.variables[type+'B_lat'][:].astype(float)
            self.lon1 = rootgrp.variables[type+'E_lon'][:].astype(float)
            self.lat1 = rootgrp.variables[type+'E_lat'][:].astype(float)
            self.len  = rootgrp.variables[type+'_len'][:].astype(float)
            self.dir  = rootgrp.variables[type+'_dir'][:].astype(float)
            self.dX   = rootgrp.variables[type+'_dX'][:].astype(float)
            self.dY   = rootgrp.variables[type+'_dY'][:].astype(float)
            t0_name   = type+'B_time'
            self.t0   = ma.array(num2date(rootgrp.variables[t0_name][:],rootgrp.variables[t0_name].units))
            t1_name   = type+'E_time'
            self.t1   = ma.array(num2date(rootgrp.variables[t1_name][:],rootgrp.variables[t1_name].units))
            if type == 'val':
                notilde = lambda x: x.replace('~','')
                id_dataset  = rootgrp.variables[type+'_id'][:]
                self.id     = list(map(notilde, chartostring(id_dataset)))
                net_dataset = rootgrp.variables[type+'_network'][:]
                self.net     = list(map(notilde, chartostring(net_dataset)))
                src_dataset = rootgrp.variables[type+'_source'][:]
                self.src     = list(map(notilde, chartostring(src_dataset)))
                # placeholder for the uncertainites and flags
                self.sX  = ma.array(np.zeros_like(self.lon0),mask=np.ones_like(self.lon0).astype('bool'))
                self.sY  = ma.array(np.zeros_like(self.lon0),mask=np.ones_like(self.lon0).astype('bool'))
                self.cXY = ma.array(np.zeros_like(self.lon0),mask=np.ones_like(self.lon0).astype('bool'))
                self.flg = ma.array(np.zeros_like(self.lon0).astype('int'),mask=np.ones_like(self.lon0).astype('bool'))
                #placeholder for relative angle
                self.relangle = ma.array(np.zeros_like(self.lon0),mask=np.ones_like(self.lon0).astype('bool'))
            else:
                self.id    = ma.array(['ID']*self.t1.shape[0])
                self.net   = ma.array(['NET']*self.t1.shape[0])
                self.src   = ma.array(['SRC']*self.t1.shape[0])
                # read the product uncertainties and flags
                self.sX  = rootgrp.variables['prd_sX'][:].astype(float)
                self.sY  = rootgrp.variables['prd_sY'][:].astype(float)
                self.cXY = rootgrp.variables['prd_cXY'][:].astype(float)
                self.flg = rootgrp.variables['prd_flag'][:].astype(float)
                # read the relative angle of orbits (if S2S)
                try:
                    self.relangle = rootgrp.variables['prd_orbit_angles'][:].astype(float)
                except KeyError:
                    self.relangle = ma.array(np.zeros_like(self.lon0),mask=np.ones_like(self.lon0).astype('bool'))


        except KeyError as e:
            raise Exception("Cannot find dataset %s in matchup file" % e)

    def _compute_extra_sid_estimates_fields(self):
        self.dt = (self.t1 - self.t0)
        dtmask = self.dt.mask
        dth = [dt.total_seconds()/(60.*60.) for dt in self.dt]
        self.dt = np.ma.array(dth,mask=dtmask)

class sid_matchup:
    def __init__(self):
        self.val  = sid_estimates()
        self.prd  = sid_estimates()
        self.coll = sid_collocation()
        self.mask = ma.nomask
        self.area = 'unknown'

    def __str__(self):
        return "VAL: {%s}, PRD: {%s}, COLL: {%s}, area:%s" % \
                (self.val.__str__(), self.prd.__str__(),self.coll.__str__(),\
                 self.area)

    def __enter__(self):
        self.val.__enter__()
        self.prd.__enter__()
        self.coll._compute(self.val, self.prd)
        self._compute_extra_matchup_fields()

    def __exit__(self, exc_type, exc_value, exc_traceback):
        self.val.__exit__(exc_type, exc_value, exc_traceback)
        self.prd.__exit__(exc_type, exc_value, exc_traceback)
        pass

    def __len__(self):
        return self.get_nb_matchups()

    def get_nb_matchups(self):
        return self.val.lon0.shape[0]

    def get_nb_selected_matchups(self):
        if self.mask is np.ma.nomask:
            return len(self)
        else:
            return (~self.mask).sum()


    def count(self):
        return self.get_nb_selected_matchups()

    def resetmask(self, mask=np.ma.nomask):

        if mask is not np.ma.nomask:
            mask = np.asarray(mask).astype('bool')
            if hasattr(mask,'mask'):
                mask = mask.filled(True)
            self.mask = np.copy(mask)
        else:
            self.mask = np.ma.nomask

        self._applymask()

    def _applymask(self):
        try:
            self.val._applymask(self.mask)
            self.prd._applymask(self.mask)
            self.coll._applymask(self.mask)
        except Exception as e:
            raise Exception("Cannot apply mask to matchup fields (%s)" % e)

    def _compute_extra_matchup_fields(self):
        # determine hemispheric area coverage of the matchup pairs
        if ((self.prd.lat0 > 0).all()):
            self.area = 'nh'
        elif ((self.prd.lat0 < 0).all()):
            self.area = 'sh'
        else:
            raise ValueError("Unable to determine area of matchups (seems global)")

    def write_to_matchup_file(self, fname):

        with Dataset(fname,'w') as _:
            mdim = _.createDimension("matchup", self.count())

            self.val.write_to_matchup_file(_)
            self.prd.write_to_matchup_file(_)


    def load_from_matchup_file(self, nc_file):
        self.load_matchups_from_list((nc_file,))

    def load_matchups_from_list(self,filenames):
        # load all matchups from given list of files
        matchups = []
        try:
            for f in filenames:
                try:
                    # read each file in a new matchup object
                    m = sid_matchup()
                    with Dataset(f) as _:
                        m.val.load_from_matchup_file(_, 'val')
                        m.prd.load_from_matchup_file(_, 'prd')
                    # append 'm' to list of matchups
                    matchups.append(m)
                except Exception as e:
                    raise(e)
                    print("Something wrong with reading file {}... skip it".format(os.path.basename(f)))
                    pass

        except Exception as e:
            raise Exception("This did not work (%s)" % e)

        # now join all the matchups to self in one go (to avoid multiple 'concatenate')
        self.val.join(list(map(attrgetter('val'),matchups)))
        self.prd.join(list(map(attrgetter('prd'),matchups)))

        # validate 'self' and compute extra matchup fields
        with self:
            pass

    def load_matchups_from_list_file(self, list_file):
        # get the list of all interesting matchup files
        try:
            filenames = read_list_of_matchup_files(list_file)
        except Exception as e:
            raise ValueError("Reading list of matchup filenames did not work (%s)" % e)
        # load them in the self object
        self.load_matchups_from_list(filenames)

    def apply_selectors(self, S, oper = 'and'):
        # select the binary logical operator to apply
        if oper == 'and':
            func = sp.ndarray.__mul__
        elif oper == 'or':
            func = sp.ndarray.__add__
        else:
            raise ValueError("Only know about oper='and', oper='or' operators")
        if self.mask is ma.nomask:
            start_select = sp.array([True] * self.get_nb_matchups())
            if oper != 'and':
                print("WARNING: The first selection on a matchup object must be oper='and'")
                oper = 'and'
                func = sp.ndarray.__mul__
        else:
            start_select = ~self.mask
        try:
            #print "Apply selector: %s" % S
            #print "Onto: %s" % msel.access_inst_field(self, S.f).data
            #print "With operator %s" % oper
            select = func(start_select,sp.array(S.select(self)))
            self.mask = ~select
            self._applymask()
        except Exception as e:
            raise Exception ("ERROR %s" % e)

    def coarse_qc_matchups(self,lim):
        """ Mask out matchups where the mismatch in dX or dY is too large,
               possibly indicating an issue with the buoys (e.g. itp79 in Oct 2014)
        """
        mismatch_x_d = self.prd.dX.filled(0) - self.val.dX.filled(0)
        mismatch_y_d = self.prd.dY.filled(0) - self.val.dY.filled(0)
        qc_mask = ((abs(mismatch_x_d) > lim) + (abs(mismatch_y_d) > lim))
        if self.mask is ma.nomask:
            new_mask = qc_mask
        else:
            new_mask = qc_mask + self.mask
        self.mask = new_mask
        self._applymask()

if __name__ == '__main__':
    m = sid_matchup()
    try:
        m.load_matchups_from_list_file('./listMatchups.txt')

        #sel1 = msel.exact_label_selector('NP-38','val.id')
        #m.apply_selectors(sel1)
        #sel2 = msel.exact_label_selector('itp43','val.id')
        #m.apply_selectors(sel2, oper='or')
        #print("Selected %s matchups (out of %s)" % (m.get_nb_selected_matchups(), m.get_nb_matchups()))
        #m.resetmask()
        #sel1 = msel.exact_label_selector('NP-38','val.id')
        #m.apply_selectors(sel1)
        #print("Selected %s matchups (out of %s)" % (m.get_nb_selected_matchups(), m.get_nb_matchups()))
    except Exception as e:
        print("ERROR: %s" % e)
        sys.exit(2)

    #print "Matchups are %s" % m
