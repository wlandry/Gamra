#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default


    restrict_msg="Checking for restrict flag "
    restrict_fragment="int main() { double * restrict x;}\n"
    restrict_flags=['--restrict','-restrict','-Drestrict=__restrict__','-Drestrict=__restrict','-Drestrict=']
    if conf.options.ftensor_restrict_flag:
        restrict_flags=[restrict_flag]

    for flag in restrict_flags:
        try:
            conf.check_cxx(msg=restrict_msg+flag, fragment=restrict_fragment,
                           cxxflags=flag, uselib_store='restrict')
        except conf.errors.ConfigurationError:
            continue
        else:
            found_restrict=True
            break

    # Find FTensor
    if conf.options.ftensor_dir:
        if not conf.options.ftensor_incdir:
            conf.options.ftensor_incdir=conf.options.ftensor_dir + "/include"

    if conf.options.ftensor_incdir:
        ftensor_incdir=[conf.options.ftensor_incdir]
    else:
        ftensor_incdir=[]

    conf.check_cxx(msg="Checking for FTensor",
                   header_name='FTensor.hpp',
                   includes=ftensor_incdir,
                   cxxflags=[flag],
                   uselib_store='FTensor')

def options(opt):
    ftensor=opt.add_option_group('FTensor Options')
    ftensor.add_option('--ftensor-dir',
                       help='Base directory where FTensor is installed')
    ftensor.add_option('--ftensor-incdir',
                       help='Directory where FTensor include files are installed')
    ftensor.add_option('--ftensor-restrict-flag',
                        help='Option to enable "restrict" in C++')
