import os
import errno

infile = 'cycle_job.sh'
outfile = 'cycle_job_tar_stgin.sh'

nnodes = 72720

wkdir = os.path.dirname(os.path.abspath(__file__))

fo = open(outfile, 'w')

with open(infile) as f:
    for line in f:
        if line[0:13] == "#PJM --stgin " or line[0:17] == "#PJM --stgin-dir ":
            res = line.find('rank=')
            reslist = line[res+5:].rstrip().rstrip("\"").split(' ')
            rank = reslist[0]
            path1 = reslist[1]
            tmppath2 = reslist[2]
            frank, path2 = tmppath2.split(':')
            dirname1 = os.path.dirname(path1)
            basename1 = os.path.basename(path1)
            dirname2 = os.path.dirname(path2)
            basename2 = os.path.basename(path2)
            if basename2 == '':
                basename2 = '.'
            if basename1 == '*':
                path1 = dirname1
                basename1 = os.path.basename(path1)
                dirname1 = os.path.dirname(path1)
                path2 = dirname2
                basename2 = os.path.basename(path2)
                dirname2 = os.path.dirname(path2)

            if rank == '*' and frank == '%r':
                print path1, dirname1, basename1, path2, dirname2, basename2

                for n in xrange(nnodes):
                    ipath = '{0:s}/input/{1:06d}/{2:s}'.format(wkdir, n, dirname2)
#                    print ipath
                    try:
                        os.makedirs(ipath)
                    except OSError as exc:  # Python >2.5
                        if exc.errno == errno.EEXIST and os.path.isdir(ipath):
                            pass
                        else:
                            raise
                    os.chdir(ipath)
                    os.symlink(path1, basename2)

            elif rank == frank:
                print path1, dirname1, basename1, path2, dirname2, basename2

                ipath = '{0:s}/input/{1:06d}/{2:s}'.format(wkdir, int(rank), dirname2)
#                print ipath
                try:
                    os.makedirs(ipath)
                except OSError as exc:  # Python >2.5
                    if exc.errno == errno.EEXIST and os.path.isdir(ipath):
                        pass
                    else:
                        raise
                os.chdir(ipath)
                os.symlink(path1, basename2)

            else:
                raise

        else:
            fo.write(line)

fo.close()

#for n in xrange(nnodes):
#    os.chdir('{0:s}/input/{1:06d}'.format(wkdir, n))
#    os.system("tar chvf ../{0:06d}.tar *".format(n))
