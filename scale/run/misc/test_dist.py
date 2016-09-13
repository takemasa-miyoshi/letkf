from distribute import *

nnodes = 2
ppn = 12
member = 23
min_np = 4
max_np = 4

#nnodes = 10
#ppn = 4
#member = 8
#min_np = 1
#max_np = 2

totalnp = ppn * nnodes

mem_nodes, mem_np, repeat_mems, parallel_mems = set_mem_np(nnodes, ppn, member, min_np, max_np)

print 'nnodes        =', nnodes
print 'ppn           =', ppn
print 'member        =', member
print 'mem_nodes     =', mem_nodes
print 'mem_np        =', mem_np
print 'repeat_mems   =', repeat_mems
print 'parallel_mems =', parallel_mems

procs, mem2node, mem2proc, proc2mem = set_mem_node_proc(nnodes, ppn, member, mem_nodes, mem_np)

print
print procs
print
i = 0
for i in xrange(member):
    print i, mem2node[i], mem2proc[i]

print
i = 0
for ip in proc2mem:
    print i, procs[i], ip
    i += 1
