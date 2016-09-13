from distribute import *

nnodes = 3636
ppn = 1
member = 100
min_np = 36
max_np = min_np

totalnp = ppn * nnodes

mem_nodes, mem_np, repeat_mems, parallel_mems = set_mem_np(nnodes, ppn, member+1, min_np, max_np)

print 'nnodes        =', nnodes
print 'ppn           =', ppn
print 'member        =', member
print 'mem_nodes     =', mem_nodes
print 'mem_np        =', mem_np
print 'repeat_mems   =', repeat_mems
print 'parallel_mems =', parallel_mems

procs, mem2node, mem2proc, proc2mem = set_mem_node_proc(nnodes, ppn, member+1, mem_nodes, mem_np)

for m in xrange(member):
    for ip, p in enumerate(mem2proc[m]):
        print "#PJM --stgin \"rank={:d} ../obsope/hist.{:04d}.pe{:06d}.nc {:d}:./\"".format(p, m+1, ip, p)

for m in xrange(member):
    for ip, p in enumerate(mem2proc[m]):
        print "#PJM --stgout \"rank={:d} {:d}:./obsda.{:04d}.{:06d}.dat ./\"".format(p, p, m+1, ip)


