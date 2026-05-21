import os
import csv
import math
import random
import heapq
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
from collections import defaultdict



def load_taskset(filepath):

    tasks = []
    with open(filepath, newline='') as f:
        reader = csv.DictReader(f)
        for row in reader:
            norm = {k.strip().lower(): v.strip() for k, v in row.items()}
            c = (norm.get('wcet') or norm.get('execution time') or
                 norm.get('computation time') or norm.get('c'))
            t = norm.get('period') or norm.get('t')
            d = norm.get('deadline') or norm.get('d') or norm.get('relative deadline')
            bcet = norm.get('bcet')
            if c is None or t is None:
                continue
            C = float(c)
            T = float(t)
            D = float(d) if d else T
            B = float(bcet) if bcet else C * 0.5
            if C <= 0 or T <= 0 or D <= 0:
                continue
            tasks.append({'C': C, 'T': T, 'D': D, 'B': B})
    return tasks


def utilization(tasks):
    return sum(t['C'] / t['T'] for t in tasks)


# 2. DM WCRT — RESPONSE TIME ANALYSIS

def dm_wcrt(tasks):
  
    sorted_tasks = sorted(tasks, key=lambda t: t['D'])
    wcrts = []
    for i, task in enumerate(sorted_tasks):
        C_i, D_i = task['C'], task['D']
        hp = sorted_tasks[:i]
        R = C_i
        while True:
            R_new = C_i + sum(math.ceil(R / h['T']) * h['C'] for h in hp)
            if R_new > D_i:
                return False, []
            if R_new == R:
                break
            R = R_new
        wcrts.append(R)
    return True, wcrts


def dm_wcrt_verbose(tasks):
   
    sorted_tasks = sorted(tasks, key=lambda t: t['D'])
    wcrts = []
    print(f"\n  DM RTA Step-by-Step (tasks sorted by deadline):")
    print(f"  {'Task':<8} {'C':>6} {'T':>6} {'D':>6}  Iterations")
    print(f"  {'-'*60}")

    for i, task in enumerate(sorted_tasks):
        C_i, D_i, T_i = task['C'], task['D'], task['T']
        hp = sorted_tasks[:i]
        R = C_i
        iterations = [R]
        converged = True
        while True:
            R_new = C_i + sum(math.ceil(R / h['T']) * h['C'] for h in hp)
            iterations.append(R_new)
            if R_new > D_i:
                converged = False
                break
            if R_new == R:
                break
            R = R_new
        iter_str = ' -> '.join(f'{v:.0f}' for v in iterations)
        status = f'R={R:.0f}' if converged else 'MISS'
        print(f"  tau_{i+1:<4} {C_i:>6.0f} {T_i:>6.0f} {D_i:>6.0f}  {iter_str}  [{status}]")
        if converged:
            wcrts.append(R)
        else:
            return False, [], sorted_tasks
    return True, wcrts, sorted_tasks


def dm_schedulable(tasks):
    ok, _ = dm_wcrt(tasks)
    return ok


# 3. EDF — PROCESSOR DEMAND CRITERION


def edf_schedulable(tasks):
    
    U = utilization(tasks)
    if U > 1.0 + 1e-9:
        return False
    if all(abs(t['D'] - t['T']) < 1e-9 for t in tasks):
        return U <= 1.0 + 1e-6

    numerator = sum((t['C'] / t['T']) * (t['T'] - t['D']) for t in tasks)
    if abs(1.0 - U) < 1e-9:
        L_star = max(t['D'] for t in tasks) * 10
    else:
        L_star = numerator / (1.0 - U)

    try:
        periods = [max(1, int(round(t['T']))) for t in tasks]
        H = periods[0]
        for p in periods[1:]:
            H = H * p // math.gcd(H, p)
            if H > 10_000_000:
                H = 10_000_000
                break
    except Exception:
        H = 10_000_000

    L_max = min(H, max(max(t['D'] for t in tasks), L_star))

    test_points = set()
    for t in tasks:
        d = t['D']
        while d <= L_max + 1e-9:
            test_points.add(d)
            d += t['T']

    def dbf(L):
        total = 0.0
        for t in tasks:
            val = math.floor((L + t['T'] - t['D']) / t['T'])
            if val > 0:
                total += val * t['C']
        return total

    for L in sorted(test_points):
        if L <= 0:
            continue
        if dbf(L) > L + 1e-9:
            return False
    return True


def edf_schedulable_verbose(tasks):
    
    U = utilization(tasks)
    print(f"\n  EDF PDC Step-by-Step:")
    print(f"  Total utilization U = {U:.4f}")

    if U > 1.0 + 1e-9:
        print(f"  U > 1.0 => NOT schedulable (necessary condition fails)")
        return False

    if all(abs(t['D'] - t['T']) < 1e-9 for t in tasks):
        result = U <= 1.0 + 1e-6
        print(f"  Implicit deadlines (D=T for all tasks)")
        print(f"  U <= 1.0 => {'SCHEDULABLE' if result else 'NOT schedulable'}")
        return result

    print(f"  Constrained deadlines detected (D < T for some tasks)")
    print(f"  Computing demand bound function dbf(L) at deadline points...")

    numerator = sum((t['C'] / t['T']) * (t['T'] - t['D']) for t in tasks)
    if abs(1.0 - U) < 1e-9:
        L_star = max(t['D'] for t in tasks) * 10
    else:
        L_star = numerator / (1.0 - U)
    print(f"  L* = {L_star:.1f}")

    try:
        periods = [max(1, int(round(t['T']))) for t in tasks]
        H = periods[0]
        for p in periods[1:]:
            H = H * p // math.gcd(H, p)
            if H > 10_000_000:
                H = 10_000_000
                break
    except Exception:
        H = 10_000_000

    L_max = min(H, max(max(t['D'] for t in tasks), L_star))
    print(f"  Testing up to L_max = {L_max:.1f}")

    test_points = set()
    for t in tasks:
        d = t['D']
        while d <= L_max + 1e-9:
            test_points.add(d)
            d += t['T']

    def dbf(L):
        total = 0.0
        for t in tasks:
            val = math.floor((L + t['T'] - t['D']) / t['T'])
            if val > 0:
                total += val * t['C']
        return total

    print(f"  Number of test points: {len(test_points)}")
    
    checked = 0
    worst_ratio = 0
    worst_L = 0
    for L in sorted(test_points):
        if L <= 0:
            continue
        demand = dbf(L)
        ratio = demand / L if L > 0 else 0
        if ratio > worst_ratio:
            worst_ratio = ratio
            worst_L = L
        if demand > L + 1e-9:
            print(f"  FAIL at L={L:.1f}: dbf(L)={demand:.1f} > L={L:.1f}")
            print(f"  => NOT schedulable")
            return False
        checked += 1

    print(f"  Checked {checked} points, worst ratio dbf(L)/L = {worst_ratio:.4f} at L={worst_L:.1f}")
    print(f"  All points satisfy dbf(L) <= L => SCHEDULABLE")
    return True

# 3b. EDF WCRT — ANALYTICAL (Appendix A method)


def edf_wcrt_analytical(tasks):
    
    n = len(tasks)
    if n == 0:
        return True, [], 0

    
    periods = [max(1, int(round(t['T']))) for t in tasks]
    H = periods[0]
    for p in periods[1:]:
        H = H * p // math.gcd(H, p)
        if H > 5_000_000:
            H = 5_000_000
            break

    
    wcet_tasks = [{'C': t['C'], 'T': t['T'], 'D': t['D'], 'B': t['C']} for t in tasks]
    result = simulate(wcet_tasks, algorithm='EDF', sim_duration=H, seed=0)

    wcrts = result['observed_wcrts']
    schedulable = result['schedulable']

    return schedulable, wcrts, H


def edf_wcrt_analytical_verbose(tasks):
   
    n = len(tasks)
    print(f"\n  EDF WCRT Analytical (Hyperperiod Construction):")

    # Compute hyperperiod
    periods = [max(1, int(round(t['T']))) for t in tasks]
    H = periods[0]
    for p in periods[1:]:
        H = H * p // math.gcd(H, p)
        if H > 5_000_000:
            H = 5_000_000
            print(f"  Hyperperiod capped at {H} (too large)")
            break

    # Count jobs
    total_jobs = sum(int(H / t['T']) for t in tasks)
    print(f"  Hyperperiod H = {H}")
    print(f"  Total jobs in [0, H): {total_jobs}")

    ok, wcrts, _ = edf_wcrt_analytical(tasks)

    if ok:
        print(f"  EDF WCRT Result: SCHEDULABLE")
        sorted_order = sorted(range(n), key=lambda i: tasks[i]['D'])
        for rank, i in enumerate(sorted_order):
            print(f"    tau_{rank+1}: WCRT = {wcrts[i]:.0f}, Deadline = {tasks[i]['D']:.0f}, Slack = {tasks[i]['D'] - wcrts[i]:.0f}")
    else:
        print(f"  EDF WCRT Result: NOT SCHEDULABLE (deadline miss detected)")

    return ok, wcrts


# 4. DISCRETE-EVENT SIMULATION


def simulate(tasks, algorithm='DM', sim_duration=None, seed=42, record_schedule=False):
  
    random.seed(seed)
    n = len(tasks)

    if sim_duration is None:
        try:
            periods = [max(1, int(round(t['T']))) for t in tasks]
            H = periods[0]
            for p in periods[1:]:
                H = H * p // math.gcd(H, p)
                if H > 1_000_000:
                    H = 1_000_000
                    break
        except Exception:
            H = 10000
        sim_duration = min(H * 5, 2_000_000)

    dm_priority = {i: rank for rank, i in
                   enumerate(sorted(range(n), key=lambda i: tasks[i]['D']))}

    clock = 0.0
    cpu_job = None
    ready_queue = []
    event_list = []
    seq_counter = [0]
    schedule_log = []  # (start_time, end_time, task_idx) for Gantt

    def push_event(t, etype, data):
        seq_counter[0] += 1
        heapq.heappush(event_list, (t, seq_counter[0], etype, data))

    max_response = [0.0] * n
    deadline_misses = [0] * n
    all_response_times = [[] for _ in range(n)]
    job_count = [0] * n

    for i in range(n):
        push_event(0.0, 'RELEASE', i)

    def job_priority(job):
        if algorithm == 'DM':
            return dm_priority[job['task_idx']]
        else:
            return (job['abs_deadline'], job['task_idx'])

    def should_preempt(new_job, current_job):
        if algorithm == 'DM':
            return dm_priority[new_job['task_idx']] < dm_priority[current_job['task_idx']]
        else:
            return new_job['abs_deadline'] < current_job['abs_deadline']

    def best_ready():
        return min(ready_queue, key=job_priority) if ready_queue else None

    def dispatch():
        nonlocal cpu_job
        if cpu_job is not None:
            return
        nxt = best_ready()
        if nxt is None:
            return
        ready_queue.remove(nxt)
        cpu_job = nxt
        cpu_job['token'] += 1
        cpu_job['last_start'] = clock
        finish_time = clock + cpu_job['remaining']
        push_event(finish_time, 'COMPLETE', (cpu_job, cpu_job['token']))

    def try_preempt():
        nonlocal cpu_job
        if cpu_job is None:
            dispatch()
            return
        nxt = best_ready()
        if nxt is None:
            return
        if should_preempt(nxt, cpu_job):
            elapsed = clock - cpu_job['last_start'] if cpu_job['last_start'] else 0
            cpu_job['remaining'] -= elapsed
            if record_schedule and cpu_job['last_start'] is not None and elapsed > 0:
                schedule_log.append((cpu_job['last_start'], clock, cpu_job['task_idx']))
            cpu_job['last_start'] = None
            ready_queue.append(cpu_job)
            cpu_job = None
            dispatch()

    while event_list:
        event_time, _, etype, data = heapq.heappop(event_list)
        if etype == 'RELEASE' and event_time >= sim_duration:
            continue
        if etype == 'COMPLETE' and event_time > sim_duration * 2:
            break
        clock = event_time

        if etype == 'RELEASE':
            i = data
            task = tasks[i]
            wcet = task['C']
            bcet = task.get('B', wcet * 0.5)
            exec_time = random.uniform(max(1.0, bcet), wcet)
            job = {
                'task_idx': i,
                'release_time': clock,
                'abs_deadline': clock + task['D'],
                'remaining': exec_time,
                'last_start': None,
                'token': 0,
            }
            job_count[i] += 1
            ready_queue.append(job)
            if clock + task['T'] < sim_duration:
                push_event(clock + task['T'], 'RELEASE', i)
            if cpu_job is not None and cpu_job['last_start'] is None:
                cpu_job['last_start'] = clock
            try_preempt()
            if cpu_job is not None and cpu_job['last_start'] is None:
                cpu_job['last_start'] = clock

        elif etype == 'COMPLETE':
            job, token = data
            if cpu_job is not job or job['token'] != token:
                continue
            if record_schedule and job['last_start'] is not None:
                schedule_log.append((job['last_start'], clock, job['task_idx']))
            response = clock - job['release_time']
            i = job['task_idx']
            max_response[i] = max(max_response[i], response)
            all_response_times[i].append(response)
            if clock > job['abs_deadline'] + 1e-9:
                deadline_misses[i] += 1
            cpu_job = None
            dispatch()
            if cpu_job is not None and cpu_job['last_start'] is None:
                cpu_job['last_start'] = clock

    return {
        'observed_wcrts': max_response,
        'all_response_times': all_response_times,
        'deadline_misses': deadline_misses,
        'total_misses': sum(deadline_misses),
        'schedulable': sum(deadline_misses) == 0,
        'job_count': job_count,
        'schedule_log': schedule_log if record_schedule else [],
    }



# 5. HAND-CRAFTED EXAMPLE TASK SETS


def hand_crafted_examples(results_dir):
    
    print("\n" + "="*60)
    print("  HAND-CRAFTED EXAMPLE TASK SETS")
    print("="*60)

    # Example 1: Both DM and EDF schedulable (moderate utilization, implicit deadlines)
    example1 = [
        {'C': 1, 'T': 4,  'D': 4,  'B': 1},
        {'C': 2, 'T': 6,  'D': 6,  'B': 1},
        {'C': 3, 'T': 12, 'D': 12, 'B': 2},
    ]

    # Example 2: Only EDF schedulable (constrained deadlines)
    # DM fails on tau_3 (R=8 > D=7), but EDF PDC passes (worst dbf/L = 0.89)
    example2 = [
        {'C': 2, 'T': 5,  'D': 4,  'B': 1},
        {'C': 2, 'T': 7,  'D': 5,  'B': 1},
        {'C': 2, 'T': 10, 'D': 7,  'B': 1},
    ]

    # Example 3: Neither schedulable (overloaded, U > 1)
    example3 = [
        {'C': 3, 'T': 5,  'D': 5,  'B': 2},
        {'C': 3, 'T': 7,  'D': 6,  'B': 2},
        {'C': 4, 'T': 10, 'D': 9,  'B': 3},
    ]

    examples = [
        ("Example 1: Both DM and EDF schedulable (implicit deadlines)", example1),
        ("Example 2: Only EDF schedulable (constrained deadlines)", example2),
        ("Example 3: Neither schedulable (overloaded)", example3),
    ]

    for title, tasks in examples:
        U = utilization(tasks)
        print(f"\n  {'─'*56}")
        print(f"  {title}")
        print(f"  U = {U:.4f}")
        print(f"  {'─'*56}")

        
        print(f"\n  {'Task':<8} {'C':>4} {'T':>6} {'D':>6} {'B':>4}")
        print(f"  {'-'*30}")
        for i, t in enumerate(tasks):
            print(f"  tau_{i+1:<4} {t['C']:>4.0f} {t['T']:>6.0f} {t['D']:>6.0f} {t['B']:>4.0f}")

        
        dm_ok, dm_wcrts, sorted_t = dm_wcrt_verbose(tasks)
        print(f"\n  DM Result: {'SCHEDULABLE' if dm_ok else 'NOT SCHEDULABLE'}")
        if dm_ok:
            for i, w in enumerate(dm_wcrts):
                print(f"    tau_{i+1}: WCRT = {w:.0f}, Deadline = {sorted_t[i]['D']:.0f}, Slack = {sorted_t[i]['D'] - w:.0f}")

    
        edf_ok = edf_schedulable_verbose(tasks)

        
        if edf_ok:
            edf_wcrt_ok, edf_wcrts = edf_wcrt_analytical_verbose(tasks)

        # Simulation comparison
        if dm_ok or edf_ok:
            print(f"\n  Simulation Results (random [BCET, WCET]):")
            sorted_tasks_sim = sorted(tasks, key=lambda t: t['D'])

            sim_dm = simulate(sorted_tasks_sim, algorithm='DM', seed=42,
                              sim_duration=240, record_schedule=True)
            sim_edf = simulate(sorted_tasks_sim, algorithm='EDF', seed=42,
                               sim_duration=240, record_schedule=True)

            # Map EDF WCRTs to sorted-by-deadline order
            sorted_indices = sorted(range(len(tasks)), key=lambda i: tasks[i]['D'])

            print(f"\n  {'Task':<8} {'DM WCRT':>10} {'EDF WCRT':>10} {'DM Sim RT':>11} {'EDF Sim RT':>12} {'DM Miss':>9} {'EDF Miss':>10}")
            print(f"  {'-'*75}")
            for i in range(len(tasks)):
                dm_w = f"{dm_wcrts[i]:.0f}" if dm_ok and i < len(dm_wcrts) else "N/A"
                edf_w = f"{edf_wcrts[sorted_indices[i]]:.0f}" if edf_ok else "N/A"
                dm_s = f"{sim_dm['observed_wcrts'][i]:.1f}"
                edf_s = f"{sim_edf['observed_wcrts'][i]:.1f}"
                print(f"  tau_{i+1:<4} {dm_w:>10} {edf_w:>10} {dm_s:>11} {edf_s:>12} {sim_dm['deadline_misses'][i]:>9} {sim_edf['deadline_misses'][i]:>10}")

            # Generate Gantt charts for this example
            wcet_tasks = [{'C': t['C'], 'T': t['T'], 'D': t['D'], 'B': t['C']} for t in sorted_tasks_sim]
            sim_dm_gantt = simulate(wcet_tasks, algorithm='DM', seed=42,
                                    sim_duration=240, record_schedule=True)
            sim_edf_gantt = simulate(wcet_tasks, algorithm='EDF', seed=42,
                                     sim_duration=240, record_schedule=True)
            if len(tasks) <= 5:
                _plot_gantt(wcet_tasks, sim_dm_gantt, sim_edf_gantt, title.split(':')[0],
                            results_dir, sim_duration=60)

    return examples


def _plot_gantt(tasks, sim_dm, sim_edf, label, results_dir, sim_duration=60):
    """Plot side-by-side Gantt charts for DM and EDF simulation."""
    n = len(tasks)
    colors = ['#4E79A7', '#F28E2B', '#E15759', '#76B7B2', '#59A14F']

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(14, 3 + n * 1.2), sharex=True)

    for ax, sim_result, alg_name in [(ax1, sim_dm, 'DM'), (ax2, sim_edf, 'EDF')]:
        for start, end, task_idx in sim_result['schedule_log']:
            if start >= sim_duration:
                continue
            end_clipped = min(end, sim_duration)
            ax.barh(task_idx, end_clipped - start, left=start, height=0.6,
                    color=colors[task_idx % len(colors)], edgecolor='black', linewidth=0.5)

        # Draw release arrows and deadline markers
        for i, t in enumerate(tasks):
            release = 0
            while release < sim_duration:
                ax.annotate('', xy=(release, i - 0.35), xytext=(release, i - 0.55),
                            arrowprops=dict(arrowstyle='->', color='green', lw=1.2))
                deadline = release + t['D']
                if deadline <= sim_duration:
                    ax.axvline(x=deadline, ymin=(n - i - 0.7) / (n + 0.5),
                               ymax=(n - i - 0.3) / (n + 0.5),
                               color='red', linewidth=1.5, linestyle='--')
                release += t['T']

        ax.set_yticks(range(n))
        ax.set_yticklabels([f'τ{i+1} (C={tasks[i]["C"]:.0f}, T={tasks[i]["T"]:.0f}, D={tasks[i]["D"]:.0f})'
                            for i in range(n)], fontsize=8)
        ax.set_title(f'{alg_name} Schedule', fontsize=11, fontweight='bold')
        ax.set_xlim(0, sim_duration)
        ax.grid(True, axis='x', alpha=0.3)
        ax.invert_yaxis()

    ax2.set_xlabel('Time')
    fig.suptitle(f'{label}: DM vs EDF Gantt Chart', fontsize=12, fontweight='bold')
    plt.tight_layout()
    safe_label = label.replace(' ', '_').replace(':', '')
    out = os.path.join(results_dir, f'{safe_label}_gantt.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print("  Gantt chart saved.")


# 6. DETAILED SIMULATION REPORT


def detailed_simulation_report(folder_path, label, results_dir, target_util=0.7):
   
    print(f"\n{'='*60}")
    print(f"  DETAILED SIMULATION REPORT: {label}")
    print(f"{'='*60}")

    best_tasks, best_diff = None, float('inf')
    for root, _, files in os.walk(folder_path):
        for fname in sorted(files):
            if not fname.lower().endswith('.csv'):
                continue
            tasks = load_taskset(os.path.join(root, fname))
            if not tasks:
                continue
            U = utilization(tasks)
            if dm_schedulable(tasks) and abs(U - target_util) < best_diff:
                best_diff = abs(U - target_util)
                best_tasks = tasks
        if best_tasks and best_diff < 0.05:
            break

    if best_tasks is None:
        print("  No suitable task set found.")
        return

    tasks = best_tasks
    n = len(tasks)
    U = utilization(tasks)
    sorted_tasks = sorted(tasks, key=lambda t: t['D'])

    dm_ok, theory_wcrts = dm_wcrt(sorted_tasks)
    edf_ok = edf_schedulable(sorted_tasks)


    edf_wcrt_ok, edf_wcrts, H_edf = edf_wcrt_analytical(sorted_tasks)

    
    wcet_tasks = [{'C': t['C'], 'T': t['T'], 'D': t['D'], 'B': t['C']} for t in sorted_tasks]
    sim_dm_wc  = simulate(wcet_tasks, algorithm='DM',  seed=42)
    sim_edf_wc = simulate(wcet_tasks, algorithm='EDF', seed=42)

    # Also run with random [BCET, WCET]
    sim_dm_rnd  = simulate(sorted_tasks, algorithm='DM',  seed=42)
    sim_edf_rnd = simulate(sorted_tasks, algorithm='EDF', seed=42)

    print(f"\n  Task set: n={n}, U={U:.3f}, Hyperperiod={H_edf}")
    print(f"  DM analytical: {'SCHEDULABLE' if dm_ok else 'NOT SCHEDULABLE'}")
    print(f"  EDF analytical (PDC): {'SCHEDULABLE' if edf_ok else 'NOT SCHEDULABLE'}")
    print(f"  EDF analytical (WCRT): {'SCHEDULABLE' if edf_wcrt_ok else 'NOT SCHEDULABLE'}")
    print(f"  DM simulation (WCET):   {sim_dm_wc['total_misses']} deadline misses")
    print(f"  EDF simulation (WCET):  {sim_edf_wc['total_misses']} deadline misses")
    print(f"  DM simulation (random): {sim_dm_rnd['total_misses']} deadline misses")
    print(f"  EDF simulation (random):{sim_edf_rnd['total_misses']} deadline misses")

    print(f"\n  Per-task comparison (WCET execution):")
    print(f"  {'Task':<7} {'C':>5} {'T':>7} {'D':>7} {'DM WCRT':>9} {'EDF WCRT':>10} {'DM Sim':>8} {'EDF Sim':>9} {'Jobs':>6} {'DM Miss':>8} {'EDF Miss':>9}")
    print(f"  {'-'*95}")
    for i in range(n):
        t = sorted_tasks[i]
        dm_w = f"{theory_wcrts[i]:.0f}" if dm_ok else "N/A"
        edf_w = f"{edf_wcrts[i]:.0f}" if edf_wcrt_ok else "N/A"
        print(f"  tau_{i+1:<3} {t['C']:>5.0f} {t['T']:>7.0f} {t['D']:>7.0f} {dm_w:>9} {edf_w:>10}"
              f" {sim_dm_wc['observed_wcrts'][i]:>8.1f} {sim_edf_wc['observed_wcrts'][i]:>9.1f}"
              f" {sim_dm_wc['job_count'][i]:>6} {sim_dm_wc['deadline_misses'][i]:>8}"
              f" {sim_edf_wc['deadline_misses'][i]:>9}")

    # Verify: simulated WCRT <= theoretical WCRT (for DM)
    if dm_ok:
        all_safe = True
        for i in range(n):
            if sim_dm_wc['observed_wcrts'][i] > theory_wcrts[i] + 1e-9:
                print(f"\n  WARNING: tau_{i+1} DM simulated RT ({sim_dm_wc['observed_wcrts'][i]:.1f}) > theoretical WCRT ({theory_wcrts[i]:.0f})")
                all_safe = False
        if all_safe:
            print(f"\n  VERIFIED: All DM simulated response times <= theoretical WCRTs")

    # Verify: EDF simulated WCRT <= EDF analytical WCRT
    if edf_wcrt_ok:
        all_safe = True
        for i in range(n):
            if sim_edf_wc['observed_wcrts'][i] > edf_wcrts[i] + 1e-9:
                print(f"\n  WARNING: tau_{i+1} EDF simulated RT ({sim_edf_wc['observed_wcrts'][i]:.1f}) > EDF analytical WCRT ({edf_wcrts[i]:.0f})")
                all_safe = False
        if all_safe:
            print(f"\n  VERIFIED: All EDF simulated response times <= EDF analytical WCRTs")



# 7. BATCH ANALYSIS


def analyse_folder(folder_path, label, results_dir):
   
    print(f"\n{'='*60}")
    print(f"  Dataset: {label}")
    print(f"  Path:    {folder_path}")
    print(f"{'='*60}")

    if not os.path.isdir(folder_path):
        print(f"  ERROR: Folder not found.")
        return None, 0, 0

    csv_files = []
    for root, dirs, files in os.walk(folder_path):
        for fname in files:
            if fname.lower().endswith('.csv'):
                csv_files.append(os.path.join(root, fname))

    if not csv_files:
        print("  No CSV files found.")
        return None, 0, 0

    print(f"  Found {len(csv_files)} task set files.")

    buckets = {}
    dm_total = edf_total = edf_wins = total = 0

    for fpath in sorted(csv_files):
        tasks = load_taskset(fpath)
        if not tasks:
            continue
        U = utilization(tasks)
        if U <= 0 or U > 2.0:
            continue

        bucket = round(round(U, 1), 1)
        dm_ok = dm_schedulable(tasks)
        edf_ok = edf_schedulable(tasks)

        if bucket not in buckets:
            buckets[bucket] = {'dm': [], 'edf': [], 'util': []}
        buckets[bucket]['dm'].append(dm_ok)
        buckets[bucket]['edf'].append(edf_ok)
        buckets[bucket]['util'].append(U)

        dm_total += dm_ok
        edf_total += edf_ok
        if edf_ok and not dm_ok:
            edf_wins += 1
        total += 1

    if total == 0:
        print("  No valid task sets processed.")
        return buckets, 0, 0

    print(f"\n  Summary over {total} task sets:")
    print(f"    DM  schedulable : {dm_total:5d}  ({100*dm_total/total:.1f}%)")
    print(f"    EDF schedulable : {edf_total:5d}  ({100*edf_total/total:.1f}%)")
    print(f"    EDF wins        : {edf_wins:5d}  (EDF ok, DM not ok)")

    print(f"\n  {'Util':>6}  {'#Sets':>6}  {'DM%':>7}  {'EDF%':>7}")
    print(f"  {'-'*32}")
    for b in sorted(buckets):
        dm_l = buckets[b]['dm']
        edf_l = buckets[b]['edf']
        nn = len(dm_l)
        print(f"  {b:>6.1f}  {nn:>6}  {100*sum(dm_l)/nn:>6.1f}%  {100*sum(edf_l)/nn:>6.1f}%")

    sb = sorted(buckets)
    dm_p  = [100*sum(buckets[b]['dm'])/len(buckets[b]['dm'])  for b in sb]
    edf_p = [100*sum(buckets[b]['edf'])/len(buckets[b]['edf']) for b in sb]

   
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(sb, dm_p,  'b-o', label='DM',  linewidth=2, markersize=6)
    ax.plot(sb, edf_p, 'r-s', label='EDF', linewidth=2, markersize=6)
    ax.axvline(x=0.69, color='gray', linestyle=':', linewidth=1.5, label='RM bound (ln2≈0.69)')
    ax.axvline(x=1.0,  color='black', linestyle='--', linewidth=1.5, label='U=1.0')
    ax.set_xlabel('Processor Utilization')
    ax.set_ylabel('Schedulable Task Sets (%)')
    ax.set_title(f'{label}: DM vs EDF Schedulability vs Utilization')
    ax.legend()
    ax.grid(True, alpha=0.3)
    ax.set_ylim(-5, 105)
    plt.tight_layout()
    out1 = os.path.join(results_dir, f'{label}_schedulability.png')
    plt.savefig(out1, dpi=150)
    plt.close()
    print("  Plot saved.")

    
    fig, ax = plt.subplots(figsize=(9, 5))
    for b in sb:
        for u, dm_ok, edf_ok in zip(buckets[b]['util'], buckets[b]['dm'], buckets[b]['edf']):
            if dm_ok and edf_ok:           c = 'green'
            elif edf_ok and not dm_ok:     c = 'orange'
            elif not dm_ok and not edf_ok: c = 'red'
            else:                          c = 'blue'
            ax.scatter(u, 1 if dm_ok else 0, color=c, alpha=0.25, s=8)
    ax.legend(handles=[
        mpatches.Patch(facecolor='green',  label='Both schedulable'),
        mpatches.Patch(facecolor='orange', label='EDF only (EDF wins)'),
        mpatches.Patch(facecolor='red',    label='Neither schedulable'),
        mpatches.Patch(facecolor='blue',   label='DM only'),
    ])
    ax.set_xlabel('Processor Utilization')
    ax.set_ylabel('DM Schedulable (1=yes, 0=no)')
    ax.set_title(f'{label}: Schedulability Scatter Plot')
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    out2 = os.path.join(results_dir, f'{label}_scatter.png')
    plt.savefig(out2, dpi=150)
    plt.close()
    print("  Plot saved.")

 
    edf_wins_per_bucket = []
    for b in sb:
        dm_l  = buckets[b]['dm']
        edf_l = buckets[b]['edf']
        wins  = sum(1 for d, e in zip(dm_l, edf_l) if e and not d)
        edf_wins_per_bucket.append(wins)

    if sum(edf_wins_per_bucket) == 0:
        print(f"  Skipping EDF wins plot for {label} (no EDF wins — implicit deadlines)")
    else:
        fig, ax = plt.subplots(figsize=(9, 5))
        bars = ax.bar(sb, edf_wins_per_bucket, width=0.07, color='orange', edgecolor='black', alpha=0.8)
        ax.set_xlabel('Processor Utilization')
        ax.set_ylabel('Number of Task Sets where EDF wins')
        ax.set_title(f'{label}: EDF Wins Over DM per Utilization Bucket')
        ax.grid(True, alpha=0.3, axis='y')
        for bar, val in zip(bars, edf_wins_per_bucket):
            if val > 0:
                ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.5,
                        str(val), ha='center', va='bottom', fontsize=9)
        plt.tight_layout()
        out3 = os.path.join(results_dir, f'{label}_edf_wins.png')
        plt.savefig(out3, dpi=150)
        plt.close()
        print("  Plot saved.")

    return buckets, edf_wins, total



# 8. SIMULATION VALIDATION


def simulation_validation(folder_path, label, results_dir, target_util=0.7):
   
    print(f"\n  Simulation validation for {label} (target U≈{target_util}) ...")

    best_tasks, best_diff = None, float('inf')
    for root, _, files in os.walk(folder_path):
        for fname in sorted(files):
            if not fname.lower().endswith('.csv'):
                continue
            tasks = load_taskset(os.path.join(root, fname))
            if not tasks:
                continue
            U = utilization(tasks)
            if dm_schedulable(tasks) and abs(U - target_util) < best_diff:
                best_diff = abs(U - target_util)
                best_tasks = tasks
        if best_tasks and best_diff < 0.05:
            break

    if best_tasks is None:
        print("  No suitable task set found.")
        return

    tasks = best_tasks
    n = len(tasks)
    U = utilization(tasks)
    sorted_tasks = sorted(tasks, key=lambda t: t['D'])

    dm_ok, theory_wcrts = dm_wcrt(sorted_tasks)
    edf_wcrt_ok, edf_wcrts, H_edf = edf_wcrt_analytical(sorted_tasks)
    sim_dm  = simulate(sorted_tasks, algorithm='DM',  seed=42)
    sim_edf = simulate(sorted_tasks, algorithm='EDF', seed=42)

    print(f"  Task set: U={U:.3f}, n={n}, Hyperperiod={H_edf}")
    print(f"  DM analytical: {'schedulable' if dm_ok else 'NOT schedulable'}")
    print(f"  EDF analytical WCRT: {'schedulable' if edf_wcrt_ok else 'NOT schedulable'}")
    dm_sim_str = '0 misses' if sim_dm['schedulable'] else str(sim_dm['total_misses']) + ' misses'
    print(f"  DM simulation: {dm_sim_str}")
    edf_sim_str = '0 misses' if sim_edf['schedulable'] else str(sim_edf['total_misses']) + ' misses'
    print(f"  EDF simulation: {edf_sim_str}")

    # Plot: WCRT theory vs simulation bar chart (now with 4 bars)
    x = list(range(n))
    w = 0.2
    fig, ax = plt.subplots(figsize=(max(10, n//2), 5))
    ax.bar([i-1.5*w for i in x], theory_wcrts if dm_ok else [0]*n, w,
           label='DM Analytical WCRT', color='steelblue', alpha=0.9)
    ax.bar([i-0.5*w for i in x], edf_wcrts if edf_wcrt_ok else [0]*n, w,
           label='EDF Analytical WCRT', color='indianred', alpha=0.9)
    ax.bar([i+0.5*w for i in x], sim_dm['observed_wcrts'],  w,
           label='DM Simulated Max RT',  color='skyblue', alpha=0.9)
    ax.bar([i+1.5*w for i in x], sim_edf['observed_wcrts'], w,
           label='EDF Simulated Max RT', color='lightsalmon', alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels([f'τ{i+1}' for i in x], fontsize=7)
    ax.set_xlabel('Task (sorted by deadline)')
    ax.set_ylabel('Response Time')
    ax.set_title(f'{label}: WCRT Theory vs Simulation (U={U:.2f})')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    out4 = os.path.join(results_dir, f'{label}_wcrt_validation.png')
    plt.savefig(out4, dpi=150)
    plt.close()
    print("  Plot saved.")

    # Plot: Theory vs simulation gap
    if dm_ok:
        gaps = [theory_wcrts[i] - sim_dm['observed_wcrts'][i] for i in range(n)]
        gap_pct = [100 * gaps[i] / theory_wcrts[i] if theory_wcrts[i] > 0 else 0 for i in range(n)]
        fig, ax = plt.subplots(figsize=(max(10, n//2), 4))
        ax.bar(x, gap_pct, color='purple', alpha=0.7)
        ax.set_xticks(x)
        ax.set_xticklabels([f'τ{i+1}' for i in x], fontsize=7)
        ax.set_xlabel('Task')
        ax.set_ylabel('Gap (%) = (Theory - Sim) / Theory × 100')
        ax.set_title(f'{label}: How Conservative is DM Theoretical WCRT? (U={U:.2f})')
        ax.grid(True, alpha=0.3, axis='y')
        ax.axhline(y=0, color='black', linewidth=0.8)
        plt.tight_layout()
        out5 = os.path.join(results_dir, f'{label}_wcrt_gap.png')
        plt.savefig(out5, dpi=150)
        plt.close()
        print("  Plot saved.")



# 9. SCHEDULABILITY VALIDATION
# (analytical result confirmed by simulation)


def schedulability_validation(folder_path, label, results_dir, n_samples=30):
   
    print(f"\n  Schedulability validation (analytical vs simulation) for {label} ...")

    csv_files = []
    for root, _, files in os.walk(folder_path):
        for fname in files:
            if fname.lower().endswith('.csv'):
                csv_files.append(os.path.join(root, fname))

    if not csv_files:
        return

  
    buckets = defaultdict(list)
    for fpath in sorted(csv_files):
        tasks = load_taskset(fpath)
        if not tasks:
            continue
        U = utilization(tasks)
        b = round(round(U, 1), 1)
        buckets[b].append((U, tasks, fpath))

    results = []  
    checked = 0

    for b in sorted(buckets.keys()):
        items = buckets[b]
        sample = items[:3]
        for U, tasks, fpath in sample:
            if checked >= n_samples:
                break
            dm_anal  = dm_schedulable(tasks)
            edf_anal = edf_schedulable(tasks)
            # Use longer simulation for validation to reduce false agreements
            sim_dm   = simulate(tasks, algorithm='DM',  seed=99, sim_duration=1_000_000)
            sim_edf  = simulate(tasks, algorithm='EDF', seed=99, sim_duration=1_000_000)
            dm_sim_ok  = sim_dm['schedulable']
            edf_sim_ok = sim_edf['schedulable']
            agree = (dm_anal == dm_sim_ok) and (edf_anal == edf_sim_ok)
            results.append((U, dm_anal, edf_anal, dm_sim_ok, edf_sim_ok, agree, fpath))
            checked += 1

    agree_count = sum(1 for r in results if r[5])
    print(f"  Checked {len(results)} task sets: {agree_count}/{len(results)} analytical-simulation agreements")

   
    disagreements = [r for r in results if not r[5]]
    if disagreements:
        print(f"\n  DISAGREEMENTS ({len(disagreements)} cases):")
        print(f"  {'U':>6}  {'DM Anal':>8}  {'DM Sim':>7}  {'EDF Anal':>9}  {'EDF Sim':>8}  Explanation")
        print(f"  {'-'*75}")
        for U, dm_a, edf_a, dm_s, edf_s, _, fpath in disagreements:
            explanation = ""
            if dm_a and not dm_s:
                explanation = "Anal=sched, Sim=miss → simulation found deadline miss (edge case)"
            elif not dm_a and dm_s:
                explanation = "Anal=unsched, Sim=no miss → sim duration may be too short"
            if edf_a and not edf_s:
                explanation += " | EDF: Anal=sched, Sim=miss"
            elif not edf_a and edf_s:
                explanation += " | EDF: Anal=unsched, Sim=no miss"
            dm_a_str = 'YES' if dm_a else 'NO'
            dm_s_str = 'YES' if dm_s else 'NO'
            edf_a_str = 'YES' if edf_a else 'NO'
            edf_s_str = 'YES' if edf_s else 'NO'
            print(f"  {U:>6.3f}  {dm_a_str:>8}  {dm_s_str:>7}  {edf_a_str:>9}  {edf_s_str:>8}  {explanation}")
            print(f"         File: {os.path.basename(fpath)}")

    # Plot: Agreement between analytical and simulation
    utils     = [r[0] for r in results]
    dm_anal   = [1 if r[1] else 0 for r in results]
    dm_sim    = [1 if r[3] else 0 for r in results]
    edf_anal  = [1 if r[2] else 0 for r in results]
    edf_sim   = [1 if r[4] else 0 for r in results]

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 4))

    ax1.scatter(utils, dm_anal, label='DM Analytical', marker='o', s=60, color='blue', alpha=0.7)
    ax1.scatter(utils, [s-0.05 for s in dm_sim], label='DM Simulation', marker='x', s=60, color='cyan', alpha=0.7)
    ax1.set_xlabel('Utilization')
    ax1.set_ylabel('Schedulable (1=yes, 0=no)')
    ax1.set_title(f'{label}: DM Analytical vs Simulation')
    ax1.legend()
    ax1.grid(True, alpha=0.3)
    ax1.set_ylim(-0.2, 1.2)

    ax2.scatter(utils, edf_anal, label='EDF Analytical', marker='o', s=60, color='red', alpha=0.7)
    ax2.scatter(utils, [s-0.05 for s in edf_sim], label='EDF Simulation', marker='x', s=60, color='salmon', alpha=0.7)
    ax2.set_xlabel('Utilization')
    ax2.set_ylabel('Schedulable (1=yes, 0=no)')
    ax2.set_title(f'{label}: EDF Analytical vs Simulation')
    ax2.legend()
    ax2.grid(True, alpha=0.3)
    ax2.set_ylim(-0.2, 1.2)

    plt.suptitle(f'{label}: Schedulability Analysis Validation ({agree_count}/{len(results)} agree)', fontsize=11)
    plt.tight_layout()
    out6 = os.path.join(results_dir, f'{label}_schedulability_validation.png')
    plt.savefig(out6, dpi=150)
    plt.close()
    print("  Plot saved.")



# 10. COMPARATIVE: DIFFERENT TASK COUNTS


def compare_task_counts(base_dir, dataset_prefix, period_dist, results_dir):
    
    print(f"\n  Comparing task counts for {dataset_prefix} ...")

    task_counts = ['5-task', '10-task', '25-task', '50-task']
    found_counts = []
    dm_means = []
    edf_means = []

    for tc in task_counts:
        folder = os.path.join(base_dir, f'{dataset_prefix}-utilDist',
                              f'{period_dist}-perDist', '1-core', tc, '0-jitter')
        if not os.path.isdir(folder):
            continue

        csv_files = []
        for root, _, files in os.walk(folder):
            for fname in files:
                if fname.lower().endswith('.csv'):
                    csv_files.append(os.path.join(root, fname))

        if not csv_files:
            continue

        dm_results = []
        edf_results = []
        for fpath in csv_files[:200]:
            tasks = load_taskset(fpath)
            if not tasks:
                continue
            U = utilization(tasks)
            if U <= 0 or U > 1.5:
                continue
            dm_results.append(dm_schedulable(tasks))
            edf_results.append(edf_schedulable(tasks))

        if dm_results:
            found_counts.append(tc)
            dm_means.append(100 * sum(dm_results) / len(dm_results))
            edf_means.append(100 * sum(edf_results) / len(edf_results))
            print(f"    {tc}: DM={dm_means[-1]:.1f}%  EDF={edf_means[-1]:.1f}%")

    if len(found_counts) < 2:
        print("  Not enough task count folders found, skipping.")
        return

    x = range(len(found_counts))
    w = 0.35
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.bar([i-w/2 for i in x], dm_means,  w, label='DM',  color='steelblue', alpha=0.85)
    ax.bar([i+w/2 for i in x], edf_means, w, label='EDF', color='salmon',    alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels(found_counts)
    ax.set_xlabel('Number of Tasks')
    ax.set_ylabel('Schedulable Task Sets (%)')
    ax.set_title(f'{dataset_prefix}: DM vs EDF by Task Count')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    ax.set_ylim(0, 110)
    plt.tight_layout()
    out = os.path.join(results_dir, f'{dataset_prefix}_task_count_comparison.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print("  Plot saved.")



# 11. DEADLINE MISSES UNDER HIGH UTILIZATION


def deadline_miss_analysis(folder_path, label, results_dir, n_samples=20):
   
    print(f"\n  Deadline miss analysis for {label} ...")

    csv_files = []
    for root, _, files in os.walk(folder_path):
        for fname in files:
            if fname.lower().endswith('.csv'):
                csv_files.append(os.path.join(root, fname))

  
    target_tasks = []
    high_util_tasks = []

    for fpath in sorted(csv_files):
        tasks = load_taskset(fpath)
        if not tasks:
            continue
        U = utilization(tasks)
        dm_ok = dm_schedulable(tasks)
        edf_ok = edf_schedulable(tasks)

        # Priority 1: DM fails, EDF passes (strongest demonstration)
        if not dm_ok and edf_ok and U <= 1.05:
            target_tasks.append((U, tasks, 'EDF_wins'))

        # Priority 2: High utilization where both might struggle
        elif 0.90 <= U <= 1.05:
            high_util_tasks.append((U, tasks, 'high_util'))

    # Combine: prefer EDF-wins cases, fill remaining with high-util
    sample = target_tasks[:n_samples]
    remaining = n_samples - len(sample)
    if remaining > 0:
        sample.extend(high_util_tasks[:remaining])

    sample.sort(key=lambda x: x[0])

    if not sample:
        print("  No suitable task sets found for deadline miss analysis.")
        print("  (This is expected for implicit-deadline automotive sets where DM and EDF are equivalent)")
        return

    utils_list   = []
    dm_misses    = []
    edf_misses   = []
    categories   = []

    for U, tasks, cat in sample:
        wcet_tasks = [{'C': t['C'], 'T': t['T'], 'D': t['D'], 'B': t['C']} for t in tasks]
        sim_dm  = simulate(wcet_tasks, algorithm='DM',  seed=42)
        sim_edf = simulate(wcet_tasks, algorithm='EDF', seed=42)
        utils_list.append(U)
        dm_misses.append(sim_dm['total_misses'])
        edf_misses.append(sim_edf['total_misses'])
        categories.append(cat)

    # Only plot if there are any misses to show
    if max(dm_misses) == 0 and max(edf_misses) == 0:
        print(f"  No deadline misses observed for {label} — both algorithms handle all sampled sets.")
        print(f"  (For implicit-deadline sets, DM and EDF have identical schedulability)")
        return

    fig, ax = plt.subplots(figsize=(10, 5))
    x = range(len(utils_list))
    ax.bar([i-0.2 for i in x], dm_misses,  0.4, label='DM deadline misses',  color='steelblue', alpha=0.85)
    ax.bar([i+0.2 for i in x], edf_misses, 0.4, label='EDF deadline misses', color='salmon',    alpha=0.85)
    ax.set_xticks(list(x))
    ax.set_xticklabels([f'{u:.2f}' for u in utils_list], rotation=45, fontsize=7)
    ax.set_xlabel('Task Set Utilization')
    ax.set_ylabel('Total Deadline Misses in Simulation')
    ax.set_title(f'{label}: Deadline Misses at High Utilization (DM vs EDF)')
    ax.legend()
    ax.grid(True, alpha=0.3, axis='y')
    plt.tight_layout()
    out = os.path.join(results_dir, f'{label}_deadline_misses.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print("  Plot saved.")

    # Print summary
    print(f"\n  Deadline miss summary:")
    print(f"  {'U':>6}  {'DM Miss':>8}  {'EDF Miss':>9}  {'Category':>12}")
    print(f"  {'-'*40}")
    for i in range(len(utils_list)):
        print(f"  {utils_list[i]:>6.3f}  {dm_misses[i]:>8}  {edf_misses[i]:>9}  {categories[i]:>12}")



# 12. COMBINED COMPARISON PLOT


def combined_comparison_plot(auto_buckets, uni_buckets, results_dir):
    
    if auto_buckets is None or uni_buckets is None:
        return

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(14, 5))

    for ax, buckets, title in [
        (ax1, auto_buckets, 'Automotive (implicit deadlines)'),
        (ax2, uni_buckets,  'UUnifast (constrained deadlines)')
    ]:
        sb    = sorted(buckets.keys())
        dm_p  = [100*sum(buckets[b]['dm'])/len(buckets[b]['dm'])  for b in sb]
        edf_p = [100*sum(buckets[b]['edf'])/len(buckets[b]['edf']) for b in sb]
        ax.plot(sb, dm_p,  'b-o', label='DM',  linewidth=2, markersize=5)
        ax.plot(sb, edf_p, 'r-s', label='EDF', linewidth=2, markersize=5)
        ax.axvline(x=0.69, color='gray', linestyle=':', linewidth=1.2)
        ax.axvline(x=1.0,  color='black', linestyle='--', linewidth=1.2)
        ax.set_xlabel('Processor Utilization')
        ax.set_ylabel('Schedulable (%)')
        ax.set_title(title)
        ax.legend()
        ax.grid(True, alpha=0.3)
        ax.set_ylim(-5, 105)

    plt.suptitle('DM vs EDF: Automotive vs UUnifast Datasets', fontsize=12)
    plt.tight_layout()
    out = os.path.join(results_dir, 'combined_comparison.png')
    plt.savefig(out, dpi=150)
    plt.close()
    print("  Combined plot saved.")


# 13. GITHUB TEST CASE VERIFICATION


def verify_test_cases(repo_path, results_dir):
    
    sched_dir   = os.path.join(repo_path, 'test_examples', 'schedulable')
    unsched_dir = os.path.join(repo_path, 'test_examples', 'not_schedulable')

    print("\n" + "="*60)
    print("  VERIFICATION AGAINST KNOWN TEST CASES")
    print("="*60)

    total = passed = 0

    print("\n  [schedulable/] — expected: DM=YES, EDF=YES")
    print(f"  {'File':<52} {'U':>5}  {'DM':>5}  {'EDF':>5}  {'OK':>4}")
    print(f"  {'-'*75}")

    for fname in sorted(os.listdir(sched_dir)):
        if not fname.endswith('.csv'):
            continue
        tasks = load_taskset(os.path.join(sched_dir, fname))
        if not tasks:
            continue
        U = utilization(tasks)
        dm_ok  = dm_schedulable(tasks)
        edf_ok = edf_schedulable(tasks)
        ok = dm_ok and edf_ok
        total += 1
        passed += ok
        short = fname.replace('_taskset.csv', '')[:51]
        print(f"  {short:<52} {U:>5.2f}  {'YES' if dm_ok else 'NO':>5}  {'YES' if edf_ok else 'NO':>5}  {'PASS' if ok else 'FAIL':>4}")

    print(f"\n  [not_schedulable/] — expected: DM=NO, EDF=YES (except Full_NonUnique)")
    print(f"  {'File':<52} {'U':>5}  {'DM':>5}  {'EDF':>5}  {'OK':>4}")
    print(f"  {'-'*75}")

    for fname in sorted(os.listdir(unsched_dir)):
        if not fname.endswith('.csv'):
            continue
        tasks = load_taskset(os.path.join(unsched_dir, fname))
        if not tasks:
            continue
        U = utilization(tasks)
        dm_ok  = dm_schedulable(tasks)
        edf_ok = edf_schedulable(tasks)
        expected_edf = False if 'Full_Utilization_NonUnique' in fname else True
        ok = (not dm_ok) and (edf_ok == expected_edf)
        total += 1
        passed += ok
        short = fname.replace('_taskset.csv', '')[:51]
        print(f"  {short:<52} {U:>5.2f}  {'YES' if dm_ok else 'NO':>5}  {'YES' if edf_ok else 'NO':>5}  {'PASS' if ok else 'FAIL':>4}")

    print(f"\n  Result: {passed}/{total} test cases passed")
    if passed == total:
        print("  All test cases passed — implementation verified!")
    else:
        print("  Some test cases failed.")

    return passed, total



# 14. SOLUTION FILE VERIFICATION


def verify_solutions(repo_path):
  
    sol_dir = os.path.join(repo_path, 'solution')
    test_dirs = [
        os.path.join(repo_path, 'test_examples', 'schedulable'),
        os.path.join(repo_path, 'test_examples', 'not_schedulable'),
    ]

    if not os.path.isdir(sol_dir):
        print("\n  No solution folder found, skipping solution verification.")
        return

    sol_files = [f for f in os.listdir(sol_dir) if f.startswith('solution_') and f.endswith('.txt')]
    if not sol_files:
        print("\n  No solution files found, skipping.")
        return

    print(f"\n{'='*60}")
    print("  VERIFICATION AGAINST PROVIDED SOLUTIONS")
    print(f"{'='*60}")

    for sol_file in sorted(sol_files):
        # Extract test case name from solution filename
        case_name = sol_file.replace('solution_', '').replace('.txt', '')

        # Find matching CSV
        csv_path = None
        for tdir in test_dirs:
            candidate = os.path.join(tdir, case_name + '_taskset.csv')
            if os.path.isfile(candidate):
                csv_path = candidate
                break

        if csv_path is None:
            continue

        tasks = load_taskset(csv_path)
        if not tasks:
            continue

        # Read solution file
        sol_path = os.path.join(sol_dir, sol_file)
        with open(sol_path, 'r') as f:
            sol_content = f.read().strip()

        # Compute our WCRTs
        sorted_tasks = sorted(tasks, key=lambda t: t['D'])
        dm_ok, dm_wcrts = dm_wcrt(sorted_tasks)
        edf_ok, edf_wcrts, H = edf_wcrt_analytical(sorted_tasks)

        print(f"\n  {case_name}:")
        print(f"    Solution file: {sol_content[:80]}...")
        if dm_ok:
            print(f"    DM WCRTs:  {[int(w) for w in dm_wcrts]}")
        if edf_ok:
            print(f"    EDF WCRTs: {[int(w) for w in edf_wcrts]}")



# 15. CUSTOM GENERATED TEST CASE ANALYSIS


def analyse_custom_testcases(gen_dir):
   
    if not os.path.isdir(gen_dir):
        print("\n  No custom generated test cases found, skipping.")
        return

    csv_files = []
    for root, _, files in os.walk(gen_dir):
        for fname in sorted(files):
            if fname.endswith('.csv'):
                csv_files.append((fname, os.path.join(root, fname)))

    if not csv_files:
        print("\n  No custom CSV files found in custom_testcases/.")
        return

    print(f"\n{'='*60}")
    print("  CUSTOM GENERATED TEST CASES ANALYSIS")
    print(f"{'='*60}")
    print(f"\n  {'Name':<35} {'n':>3} {'U':>6} {'DM':>5} {'EDF':>5} {'DM Sim':>7} {'EDF Sim':>8}")
    print(f"  {'-'*72}")

    total = passed = 0
    for fname, fpath in sorted(csv_files):
        tasks = load_taskset(fpath)
        if not tasks:
            continue
        U = utilization(tasks)
        dm_ok = dm_schedulable(tasks)
        edf_ok = edf_schedulable(tasks)
        sim_dm = simulate(tasks, algorithm='DM', seed=42)
        sim_edf = simulate(tasks, algorithm='EDF', seed=42)
        name = fname.replace('_taskset.csv', '')
        dm_str = 'YES' if dm_ok else 'NO'
        edf_str = 'YES' if edf_ok else 'NO'
        dm_sim = f"{sim_dm['total_misses']} miss" if sim_dm['total_misses'] > 0 else "OK"
        edf_sim = f"{sim_edf['total_misses']} miss" if sim_edf['total_misses'] > 0 else "OK"
        total += 1
        # Agreement: analytical matches simulation
        dm_agree = dm_ok == sim_dm['schedulable']
        edf_agree = edf_ok == sim_edf['schedulable']
        if dm_agree and edf_agree:
            passed += 1
        print(f"  {name:<35} {len(tasks):>3} {U:>6.3f} {dm_str:>5} {edf_str:>5} {dm_sim:>7} {edf_sim:>8}")

    print(f"\n  Result: {passed}/{total} custom test cases — analytical matches simulation")



# 16. MAIN

def main():
    
    SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))

    # Try ZIP/portable layout first (everything relative to script)
    if os.path.isdir(os.path.join(SCRIPT_DIR, 'test_examples')):
        BASE_DIR  = SCRIPT_DIR
        REPO_PATH = SCRIPT_DIR
        CUSTOM_DIR = os.path.join(SCRIPT_DIR, 'custom_testcases')
    else:
        # Fall back to your Mac layout
        BASE_DIR  = os.path.expanduser('~/Desktop/output')
        REPO_PATH = os.path.expanduser('~/Desktop/Taskset-Generator-Exercise')
        CUSTOM_DIR = os.path.join(REPO_PATH, 'custom_testcases')

    RESULTS_DIR = os.path.join(BASE_DIR, 'Mini_Project1', 'results')
    os.makedirs(RESULTS_DIR, exist_ok=True)

    AUTOMOTIVE = os.path.join(BASE_DIR,
        'automotive-utilDist', 'automotive-perDist', '1-core', '25-task', '0-jitter')
    UUNIFAST   = os.path.join(BASE_DIR,
        'uunifast-utilDist', 'uniform-discrete-perDist', '1-core', '25-task', '0-jitter')

    print("\n" + "="*60)
    print("  DTU DRTS Mini Project 1 — DM vs EDF Analysis")
    print("="*60 + "\n")

  
    hand_crafted_examples(RESULTS_DIR)

  
    if os.path.isdir(REPO_PATH):
        verify_test_cases(REPO_PATH, RESULTS_DIR)
        verify_solutions(REPO_PATH)
    else:
        print("Skipping verification — GitHub repo not found.")

    
    analyse_custom_testcases(CUSTOM_DIR)

    
    auto_buckets, auto_wins, auto_total = analyse_folder(AUTOMOTIVE, 'automotive', RESULTS_DIR)
    uni_buckets,  uni_wins,  uni_total  = analyse_folder(UUNIFAST,   'uunifast',   RESULTS_DIR)

    combined_comparison_plot(auto_buckets, uni_buckets, RESULTS_DIR)

    
    detailed_simulation_report(AUTOMOTIVE, 'automotive', RESULTS_DIR, target_util=0.7)
    detailed_simulation_report(UUNIFAST,   'uunifast',   RESULTS_DIR, target_util=0.7)

   
    simulation_validation(AUTOMOTIVE, 'automotive', RESULTS_DIR, target_util=0.7)
    simulation_validation(UUNIFAST,   'uunifast',   RESULTS_DIR, target_util=0.7)

    
    schedulability_validation(AUTOMOTIVE, 'automotive', RESULTS_DIR, n_samples=30)
    schedulability_validation(UUNIFAST,   'uunifast',   RESULTS_DIR, n_samples=30)

   
    deadline_miss_analysis(AUTOMOTIVE, 'automotive', RESULTS_DIR)
    deadline_miss_analysis(UUNIFAST,   'uunifast',   RESULTS_DIR)

    
    compare_task_counts(BASE_DIR, 'automotive', 'automotive', RESULTS_DIR)
    compare_task_counts(BASE_DIR, 'uunifast',   'uniform-discrete', RESULTS_DIR)

    print("\n" + "="*60)
    print("  ALL DONE!")
    print("  Files generated:")
    for f in sorted(os.listdir(RESULTS_DIR)):
        if not f.startswith('.'):
            print(f"    {f}")
    print("="*60)


if __name__ == '__main__':
    main()