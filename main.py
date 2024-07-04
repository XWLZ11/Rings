import numpy as np
import MDAnalysis as mda
from HBond import HBond
import os, time, json, sys
import networkx as nx
from itertools import groupby

from shared import rings as rg

import multiprocessing, queue


def main():
    dir_source ="/home/tsingularity/Desktop/MD/SCAN/cc.ct/CC/"
    filename = "traj.lammpstrj"

    #----------------------必要的基本设置--------------------------
    parameter = {"start": 0,
                 "stop": 1,
                 "step": 1,
                 "size_limit": 9,
                 "definition": 2,
                 "num_O_atoms": 129,   #体系氧原子数量
                 }


    #------------------------指定并行数----------------------------
    if len(sys.argv) < 2:
        print("Usage: python script.py <integer>")
        sys.exit(1)
    try:
        parameter["np"] = int(sys.argv[1])
    except ValueError:
        print("The argument must be an integer.")
        sys.exit(1)
    st = time.time()
    rings(dir_source, filename, param=parameter)
    ed = time.time()
    print(f"Time used: {ed-st:.2f}s")

def rings(dir_source, filename, param=None):
    #-------------------------氢键定义---------------------------------
    u = mda.Universe(os.path.join(dir_source, filename), topology_format="LAMMPSDUMP", format="LAMMPSDUMP")

    hbonds_don = HBond(
        universe=u,
        donors_sel='index 126 127 128',
        hydrogens_sel='type 2',
        acceptors_sel='index 0:125',
        d_h_cutoff=1.3,
        d_a_cutoff=3.0,
        h_d_a_angle_cutoff=25,
        update_selections=False)
    hbonds_acc = HBond(
        universe=u,
        donors_sel='index 0:125',
        hydrogens_sel='type 2',
        acceptors_sel='index 126 127 128',
        d_h_cutoff=1.3,
        d_a_cutoff=3.0,
        h_d_a_angle_cutoff=25,
        update_selections=False)
    hbonds_water = HBond(
        universe=u,
        donors_sel='index 0:125',
        hydrogens_sel='type 2',
        acceptors_sel='index 0:125',
        d_h_cutoff=1.3,
        d_a_cutoff=3.5,
        h_d_a_angle_cutoff=30,
        update_selections=False)

    #------------------------环统计-----------------------------------

    hbonds_don.run(start=param["start"], stop=param["stop"], step=param["step"], verbose=True)
    hbonds_acc.run(start=param["start"], stop=param["stop"], step=param["step"], verbose=True)
    hbonds_water.run(start=param["start"], stop=param["stop"], step=param["step"], verbose=True)
    hbonds = np.vstack([hbonds_don.results['hbonds'], hbonds_acc.results['hbonds'], hbonds_water.results['hbonds']])

    result_dict = mpi_run(rings_stats, param, u, hbonds)
    dict_rings = dict(sorted(result_dict.items()))


    #-----------------------结果存为json-----------------------------
    save_json(dir_source, dict_rings)



def mpi_run(func, param, u, hbonds):
    num_frames = (param["stop"] - param["start"]) // param["step"]
    num_processes = param["np"]   # 获取可用 CPU 核心数

    # 计算每个进程需要处理的帧数
    frames_per_process = num_frames // num_processes
    remaining_frames = num_frames % num_processes

    # 用于结果收集
    manager = multiprocessing.Manager()
    result_dict = manager.dict()

    # 创建一个线程安全的队列来存储位置和步骤信息
    step_queues = [queue.Queue() for i in range(num_processes)]
    positions_queues = [queue.Queue() for i in range(num_processes)]

    # 将帧的索引放入队列中
    start_idx = param["start"]

    for i in range(num_processes):
        end_idx = start_idx + frames_per_process
        if i < remaining_frames:
            end_idx += 1  # 将多余的帧分配给前面的进程
        # 将位置和步骤信息放入队列中
        for step, ts in zip(range(start_idx, end_idx, param["step"]), u.trajectory[start_idx:end_idx:param["step"]]):
            positions_queues[i].put(u.select_atoms("all").positions)
            step_queues[i].put(step)
        start_idx = end_idx

    processes = []
    # 创建多个进程来处理数据
    for _ in range(num_processes):
        process = multiprocessing.Process(target=func, 
                                          args=(step_queues[_], 
                                                positions_queues[_], 
                                                hbonds, u, 
                                                result_dict, 
                                                param["size_limit"], 
                                                param["definition"], 
                                                param["num_O_atoms"]))
        process.start()
        processes.append(process)

    # 等待所有进程完成
    for process in processes:
        process.join()

    return result_dict


def rings_stats(step_queue, positions_queue, 
                hbonds, u, result_dict, 
                size_limit, definition, 
                num_O_atoms) -> None:

    while True:
        try:
            step = step_queue.get_nowait()
            positions = positions_queue.get_nowait()
        except queue.Empty:
            break
        else:
            frame = hbonds[:, 0] == step

            bonds = hbonds[frame, 1:4:2].astype(np.uint16)

            rings_single = networkx_ring(bonds, positions, u.dimensions[:3], size_limit=size_limit, definition=definition, num_O_atoms=num_O_atoms)
            result_dict[step] = rings_single

            print(f"Frame:    {step}\n----------------------------\n")
            for k, v in rings_single.items():
                print(f"size:    {k}, num:    {len(v)}\n")
                print(v)
                print("\n")


def networkx_ring(bonds, pos, celldm, 
                  size_limit=6, definition=2, 
                  num_O_atoms=128) -> dict:
    if size_limit <= 11:    # 选择C++或Python搜寻环函数
        Cpp_findrings = True  
    else:
        Cpp_findrings = False
    Cpp_validrings = False  # 选择C++或Python验证环函数，两者相差不多

    if Cpp_findrings:
        # 小环快(小于等于11元环时)
        # Cpp深度优先算法（找环非常全，但很多是大环包含小环，这些是被我们排除的）
        pre_rings = rg.findrings(bonds, num_O_atoms, size_limit)
    else: # 图论约翰逊算法（大于11元环时快）
        G = nx.Graph()
        for bond in bonds:
            donor_oxygen = bond[0]
            acceptor_oxygen = bond[1]
            G.add_edge(donor_oxygen, acceptor_oxygen)
        pre_rings = sorted(nx.simple_cycles(G, length_bound=size_limit))




    sum_rings = {size: np.empty((0, size), dtype=np.uint16) for size in range(3, size_limit + 1)}
    groups = groupby(sorted(pre_rings, key=len), len)
    unsum_cycles = {k: np.array(list(v), dtype=np.uint16) for k, v in groups}   #未考虑周期性的环/链

    for k, v in unsum_cycles.items():   #排除链
        for ring in v:
            if sumUp(pos, ring, celldm): # 周期性边界条件
                sum_rings[k] = np.append(sum_rings[k], np.atleast_2d(ring), axis=0)


    # 验证环
    if Cpp_validrings:  # validrings(sum_rings, definition) 纯C++实现的环验证函数，定义见shared/validring.h
        rings_single = rg.validrings(sum_rings, definition)
    else:  # 与validing等价的python程序，性能方面相差无几，性能依赖于new_bond函数的性能，但已经用C++优化
        [rings_single, rings_single_all] = [{size: np.empty((0, size), dtype=np.uint16) for size in range(3, size_limit + 1)} for i in range(2)]
        for k, v in sum_rings.items():
            tmp = []
            for ring in v:
                if rg.new_bond(ring, rings_single, rings_single_all, definition):
                    tmp.append(ring)
            rings_single_all[k] = np.append(rings_single_all[k], v, axis=0)
            if len(tmp) == 0:
                continue
            rings_single[k] = np.append(rings_single[k], np.array(tmp), axis=0)

    return rings_single


def sumUp(pos, ring, celldm):
    ring1 = ring.copy()
    ring2 = np.roll(ring1, -1)

    ds = [pos[i1, :] - pos[i2, :] for i1, i2 in zip(ring1, ring2)]
    ds = np.array(ds)


    bool_xp = ds[:, 0] < -0.5 * celldm[0]
    bool_xn = ds[:, 0] > 0.5 * celldm[0]
    bool_yp = ds[:, 1] < -0.5 * celldm[1]
    bool_yn = ds[:, 1] > 0.5 * celldm[1]
    bool_zp = ds[:, 2] < -0.5 * celldm[2]
    bool_zn = ds[:, 2] > 0.5 * celldm[2]

    for i, mask in zip(range(3), [bool_xp, bool_yp, bool_zp]):
        ds[mask, i] += celldm[i]
    for j, mask in zip(range(3), [bool_xn, bool_yn, bool_zn]):
        ds[mask, j] -= celldm[j]

    if np.max([np.abs(np.sum(ds[:, 0])), np.abs(np.sum(ds[:, 1])), np.abs(np.sum(ds[:, 2]))]) > 1.0e-3:
        return False
    else:
        return True

def save_json(dir_source, dict_rings):
    def ndarray2list(arr):
        if isinstance(arr, np.ndarray):
            return arr.tolist()
        else:
            return arr

    json_data = json.dumps(dict_rings, default=ndarray2list, indent=0)

    with open(dir_source+'rings.json', 'w') as fp:
        fp.write(json_data)




if __name__ == "__main__":
    import warnings
    warnings.filterwarnings("ignore")
    main()
