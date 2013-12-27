#! /usr/bin/env python
# encoding: utf-8

def configure(conf):
    def get_param(varname,default):
        return getattr(Options.options,varname,'')or default

    conf.load('FTensor')
    # Find OKADA
    if conf.options.okada_dir:
        if not conf.options.okada_incdir:
            conf.options.okada_incdir=conf.options.okada_dir + "/include"
        if not conf.options.okada_libdir:
            conf.options.okada_libdir=conf.options.okada_dir + "/lib"

    if conf.options.okada_incdir:
        okada_incdir=[conf.options.okada_incdir]
    else:
        okada_incdir=[]
    if conf.options.okada_libdir:
        okada_libdir=[conf.options.okada_libdir]
    else:
        okada_libdir=[]

    if conf.options.okada_libs:
        okada_libs=conf.options.okada_libs.split()
    else:
        okada_libs=['okada']

    conf.check_cxx(msg="Checking for Okada",
                   header_name='Okada.hxx',
                   includes=okada_incdir,
                   use=['FTensor'],
                   uselib_store='okada',
                   libpath=okada_libdir,
                   rpath=okada_libdir,
                   lib=okada_libs)

def options(opt):
    opt.load('FTensor')
    okada=opt.add_option_group('Okada Options')
    okada.add_option('--okada-dir',
                   help='Base directory where okada is installed')
    okada.add_option('--okada-incdir',
                   help='Directory where okada include files are installed')
    okada.add_option('--okada-libdir',
                   help='Directory where okada library files are installed')
    okada.add_option('--okada-libs',
                   help='Names of the okada libraries without prefix or suffix\n'
                   '(e.g. "okada")')
