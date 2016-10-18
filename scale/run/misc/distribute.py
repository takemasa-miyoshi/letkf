def set_mem_np(nnodes, ppn, member, min_np=None, max_np=None):
    if nnodes >= member:
        mem_nodes = nnodes / member
        mem_np = ppn * mem_nodes
    else:
        mem_nodes = 1
        mempn = (member - 1) / nnodes + 1
        if mempn > ppn:
            mem_np = 1
        else:
            mem_np = ppn / mempn

    if min_np is not None and mem_np < min_np:
        mem_np = min_np
    if max_np is not None and mem_np > max_np:
        mem_np = max_np

    mem_nodes = (mem_np - 1) / ppn + 1
    if mem_nodes > nnodes:
        raise ValueError, '# of nodes is insufficient.'

    repeat_mems = nnodes / mem_nodes
    if mem_nodes > 1:
        parallel_mems = repeat_mems
    else:
        parallel_mems = repeat_mems * (ppn / mem_np)

    return mem_nodes, mem_np, repeat_mems, parallel_mems


def set_mem_node_proc(nnodes, ppn, member, mem_nodes, mem_np):
    procs = []
    for n in xrange(nnodes):
        procs += [n] * ppn

    if mem_nodes > 1:
        n_mem = nnodes / mem_nodes
        n_mempn = 1
    else:
        n_mem = nnodes
        n_mempn = ppn / mem_np
    nitmax = (member - 1) / (n_mem * n_mempn) + 1
    tppn = mem_np / mem_nodes
    tmod = mem_np % mem_nodes

    mem2node = [None] * member
    mem2proc = [None] * member
    proc2mem = [[None for i in range(nitmax)] for j in range(ppn * nnodes)]
    m = 0
    for it in xrange(nitmax):
        for i in xrange(n_mempn):
            n = 0
            for j in xrange(n_mem):
                if m < member:
                    mem2node[m] = []
                    mem2proc[m] = []
                    qs = 0
                    for nn in xrange(mem_nodes):
                        if nn < tmod:
                            tppnt = tppn + 1
                        else:
                            tppnt = tppn
                        for q in xrange(tppnt):
                            ip = (n+nn)*ppn + i*mem_np + q
                            mem2node[m] += [n+nn]
                            mem2proc[m] += [ip]
                            proc2mem[ip][it] = [m, qs]
                            qs += 1

                    m += 1
                    n += mem_nodes

    return procs, mem2node, mem2proc, proc2mem
