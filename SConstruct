import os, glob

cenv = Environment(ENV = {'PATH':os.environ['PATH']})
cenv.Append(CFLAGS = ['-fno-leading-underscore'])
cenv.Object('src/unixtime.c')

env = Environment(ENV = {'PATH':os.environ['PATH']})
#env.Append(FORTRANPATH = ['/group/clas/builds/centos65/include/'], LIBPATH = ['/group/clas/builds/centos65/lib/'],LIBS = ['bosio','seb','bankdefs','c_cern','clasutil'])
#env.Append(FORTRANPATH = ['/group/clas/builds/old/centos62/trunk/build/include/'], LIBPATH = ['/group/clas/builds/old/centos62/trunk/build/lib/'],LIBS = ['bosio','seb','bankdefs','c_cern','clasutil'])
#env.Append(FORTRANPATH = ['/apps/cernlib/x86_64_rhel6_4.7.2/2005/include/'], LIBPATH = ['/apps/cernlib/x86_64_rhel6_4.7.2/2005/lib/'],LIBS = ['pawlib','packlib'])
env.Append(FORTRANPATH = ['/apps/cernlib/x86_64_rhel7/2005/include/'], LIBPATH = ['/apps/cernlib/x86_64_rhel7/2005/lib/'],LIBS = ['pawlib','packlib'])
env.Append(FORTRANFLAGS = ['-O3','-m64','-fno-automatic','-ffixed-line-length-none','-fno-second-underscore'])
env.Append(FORTRANPATH = ['inc/'])

#sources = glob.glob('bos*.F')
#sources += glob.glob('h2*.F')
#sources += glob.glob('include/*.inc')
sources = ['src/unixtime.o','src/elast_gen.F']
print "sources =", sources

env.Program('elast_gen',sources)

