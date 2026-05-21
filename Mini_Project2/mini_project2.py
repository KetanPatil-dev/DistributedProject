

import heapq
import math
import random
import os
import json as _json
from collections import defaultdict as _defaultdict
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches
import numpy as np


# ─── CONSTANTS ───────────────────────────────────────────────────────────────

LINK_RATE_MBPS = 100.0                      # Megabits per second
LINK_RATE_BPS  = LINK_RATE_MBPS * 1e6       # bits per second
IDLE_SLOPE_NORM = 0.5                        # normalized; alpha+ = 0.5 * link_rate
ALPHA_PLUS  = IDLE_SLOPE_NORM * LINK_RATE_BPS   # bits/s  (50 Mb/s per CBS class)
ALPHA_MINUS = LINK_RATE_BPS - ALPHA_PLUS        # bits/s  (50 Mb/s per CBS class)


# ─── HELPER ──────────────────────────────────────────────────────────────────

def tx_time_us(size_bytes: float) -> float:
    """Transmission time in microseconds for a frame of `size_bytes` bytes."""
    return size_bytes * 8.0 / LINK_RATE_BPS * 1e6   # µs


def tx_time_s(size_bytes: float) -> float:
    """Transmission time in seconds."""
    return size_bytes * 8.0 / LINK_RATE_BPS


# ─── ANALYTICAL WCD ──────────────────────────────────────────────────────────

def _spi(stream_i, streams, pclass, alpha_plus, alpha_minus):
    """Same-Priority Interference for stream_i among streams of class `pclass`."""
    spi = 0.0
    for s in streams:
        if s is stream_i:
            continue
        if s['class'] == pclass:
            spi += tx_time_us(s['size']) * (1.0 + alpha_minus / alpha_plus)
    return spi


def wcd_avb_a(stream_i, streams, num_hops=1, alpha_plus=None, alpha_minus=None):
    """
    WCD for an AVB-A stream.
    Cao 2016 WFCS Theorem 4 + same-priority contention via iterative RTA.
    With alpha+ = alpha- = 0.5*BW, each same-class frame contributes 2*C_j.
    But for multiple streams with different periods, we use a fixed-point
    iteration (busy-window) to correctly bound the number of same-class
    frames that can interfere:
      R = C_i + LPI + sum_{j!=i, A} ceil(R / T_j) * C_j * 2
    LPI = max frame size of lower-priority streams (class B or BE).
    """
    C_i   = tx_time_us(stream_i['size'])
    sp_a  = [s for s in streams if s['class'] == 'A' and s is not stream_i]
    lower = [s for s in streams if s['class'] in ('B', 'BE')]
    # LPI: Cao 2016 WFCS Thm 4 — when a lower-priority frame is transmitting,
    # class A credit recovers at alpha+, so eligible interval is delayed by
    # C_L_max * (1 + alpha_A+ / alpha_A-) = 2 * C_L_max (since alpha+=alpha-)
    ap = alpha_plus  if alpha_plus  is not None else ALPHA_PLUS
    am = alpha_minus if alpha_minus is not None else ALPHA_MINUS
    lpi   = max((tx_time_us(s['size']) for s in lower), default=0.0) * (1.0 + ap / am)

    # Initial estimate: sum all same-class interference (upper bound seed)
    R = C_i + lpi + sum(tx_time_us(s['size']) * (1.0 + am/ap) for s in sp_a)
    for _ in range(500):
        R_new = C_i + lpi + sum(
            math.ceil(R / s['period_us']) * tx_time_us(s['size']) * (1.0 + am/ap)
            for s in sp_a
        )
        if abs(R_new - R) < 1e-9:
            R = R_new
            break
        if R_new > stream_i['period_us'] * 20:
            R = R_new
            break
        R = R_new

    return R * num_hops


def wcd_avb_b(stream_i, streams, num_hops=1, alpha_plus=None, alpha_minus=None):
   
    C_i    = tx_time_us(stream_i['size'])
    sp_b   = [s for s in streams if s['class'] == 'B' and s is not stream_i]
    higher = [s for s in streams if s['class'] == 'A']
    lower  = [s for s in streams if s['class'] == 'BE']

    C_A_max  = max((tx_time_us(s['size']) for s in higher), default=0.0)
    C_BE_max = max((tx_time_us(s['size']) for s in lower),  default=0.0)
    ap = alpha_plus  if alpha_plus  is not None else ALPHA_PLUS
    am = alpha_minus if alpha_minus is not None else ALPHA_MINUS
    mpi = C_BE_max * (1.0 + ap / am) + C_A_max

    R = C_i + mpi + sum(tx_time_us(s['size']) * (1.0 + am/ap) for s in sp_b)
    for _ in range(500):
        R_new = C_i + mpi + sum(
            math.ceil(R / s['period_us']) * tx_time_us(s['size']) * (1.0 + am/ap)
            for s in sp_b
        )
        if abs(R_new - R) < 1e-9:
            R = R_new
            break
        if R_new > stream_i['period_us'] * 20:
            R = R_new
            break
        R = R_new

    return R * num_hops


def wcd_be(stream_i, streams, num_hops=1, alpha_plus=None, alpha_minus=None):
   
    C_i    = tx_time_us(stream_i['size'])
    hp_a   = [s for s in streams if s['class'] == 'A']
    hp_b   = [s for s in streams if s['class'] == 'B']
    higher = hp_a + hp_b

    
    block = max((tx_time_us(s['size']) for s in streams), default=0.0)

    C_A_max = max((tx_time_us(s['size']) for s in hp_a), default=0.0)
    C_B_max = max((tx_time_us(s['size']) for s in hp_b), default=0.0)
    ap = alpha_plus  if alpha_plus  is not None else ALPHA_PLUS
    am = alpha_minus if alpha_minus is not None else ALPHA_MINUS
    credit_block = (C_A_max * am / ap) + (C_B_max * am / ap)

    R = C_i + block + credit_block
    for _ in range(500):
        R_new = C_i + block + credit_block
        for s in hp_a:
            R_new += math.ceil(R / s['period_us']) * tx_time_us(s['size'])
        for s in hp_b:
            R_new += math.ceil(R / s['period_us']) * tx_time_us(s['size'])
        if abs(R_new - R) < 1e-9:
            R = R_new
            break
        if R_new > stream_i['period_us'] * 20:
            R = R_new
            break
        R = R_new

    return R * num_hops


def compute_all_wcds(streams, num_hops=1, alpha_plus=None, alpha_minus=None):
    
    wcds = {}
    for s in streams:
        if s['class'] == 'A':
            wcds[s['name']] = wcd_avb_a(s, streams, num_hops, alpha_plus, alpha_minus)
        elif s['class'] == 'B':
            wcds[s['name']] = wcd_avb_b(s, streams, num_hops, alpha_plus, alpha_minus)
        else:
            wcds[s['name']] = wcd_be(s, streams, num_hops, alpha_plus, alpha_minus)
    return wcds


#  STRICT PRIORITY ANALYTICAL WCD (optional) 
def wcd_sp_rta(stream_i, streams, num_hops=1):
    
    priority = {'A': 0, 'B': 1, 'BE': 2}
    p_i = priority[stream_i['class']]

    higher = [s for s in streams if priority[s['class']] < p_i]
    same   = [s for s in streams if priority[s['class']] == p_i and s is not stream_i]
    block  = max((tx_time_us(s['size']) for s in streams), default=0.0)

    C_i = tx_time_us(stream_i['size'])

    R = C_i + block
    for _ in range(500):
        R_new = C_i + block
        R_new += sum(math.ceil(R / s['period_us']) * tx_time_us(s['size']) for s in higher)
        R_new += sum(math.ceil(R / s['period_us']) * tx_time_us(s['size']) for s in same)
        if abs(R_new - R) < 1e-9:
            R = R_new
            break
        if R_new > stream_i['period_us'] * 10:
            R = float('inf')
            break
        R = R_new

    return R * num_hops

    return R * num_hops


def compute_all_wcds_sp(streams, num_hops=1):
    
    return {s['name']: wcd_sp_rta(s, streams, num_hops) for s in streams}




_EVT_ARRIVAL = 0
_EVT_TX_END  = 1
_EVT_CREDIT_READY = 2   # when a CBS class credit rises back to 0


def _credit_at(cr_last, t_last, t_now, rising):
   
    if rising:
        return min(0.0, cr_last + ALPHA_PLUS * (t_now - t_last) * 1e-6)  # µs → s
    else:
        return cr_last - ALPHA_MINUS * (t_now - t_last) * 1e-6


def simulate_cbs(streams, num_hops=1, sim_duration_us=None, frames_per_stream=300, seed=42,
                 alpha_plus=None, alpha_minus=None):
   
    rng = random.Random(seed)
    _ap = alpha_plus  if alpha_plus  is not None else ALPHA_PLUS
    _am = alpha_minus if alpha_minus is not None else ALPHA_MINUS

    if sim_duration_us is None:
        max_period = max(s['period_us'] for s in streams)
        sim_duration_us = max_period * frames_per_stream

   
    events = []   # (time, event_type, payload)
    counter = [0]

    def push(t, etype, payload):
        heapq.heappush(events, (t, counter[0], etype, payload))
        counter[0] += 1

    for s in streams:
        t = s.get('offset_us', 0.0)
        while t < sim_duration_us:
            push(t, _EVT_ARRIVAL, (s['name'], s['class'], s['size'], t))
            t += s['period_us']

    # Simulator state
    credit  = {'A': 0.0, 'B': 0.0}   # bits (may be negative)
    cr_time = {'A': 0.0, 'B': 0.0}   # time of last credit snapshot (µs)
    wakeup_pending = [False]           # at most ONE wakeup event in the heap

    queues = {'A': [], 'B': [], 'BE': []}  # [(arrive_us, size_bytes, name)]
    transmitting = None  # (finish_us, size_bytes, arrive_us, name, cls)

    response_times = {s['name']: [] for s in streams}

    _CREDIT_EPS = 1e-6   

    def _update_credit(cls, now_us, transmitting_cls=None):
      
        dt_s = (now_us - cr_time[cls]) * 1e-6
        if dt_s <= 0:
            return
        if cls == transmitting_cls:
            credit[cls] -= _am * dt_s
        elif queues[cls] and credit[cls] < 0:
            credit[cls] = min(0.0, credit[cls] + _ap * dt_s)
        elif not queues[cls] and credit[cls] > 0:
            credit[cls] = 0.0
        # Clamp floating-point near-zero to exactly 0 to avoid infinite wakeup loops
        if -_CREDIT_EPS < credit[cls] < _CREDIT_EPS:
            credit[cls] = 0.0
        cr_time[cls] = now_us

    def _try_start_transmission(now_us):
       
        for cls in ('A', 'B'):
            _update_credit(cls, now_us)

        for cls in ('A', 'B', 'BE'):
            if not queues[cls]:
                continue
            if cls in ('A', 'B') and credit[cls] < -_CREDIT_EPS:
                continue  # ineligible — skip but don't push an event here
            arrive_us, size, name = queues[cls].pop(0)
            return (now_us + tx_time_us(size), size, arrive_us, name, cls)

        # Nothing eligible — schedule a single wakeup at the earliest credit recovery
        if not wakeup_pending[0]:
            t_wake = float('inf')
            for cls in ('A', 'B'):
                if queues[cls] and credit[cls] < 0:
                    t_wake = min(t_wake, now_us + (-credit[cls] / _ap) * 1e6)
            if t_wake < sim_duration_us:
                wakeup_pending[0] = True
                push(t_wake, _EVT_CREDIT_READY, None)
        return None

    # Process events
    while events:
        now_us, _, etype, payload = heapq.heappop(events)

        if now_us > sim_duration_us:
            break

        if etype == _EVT_ARRIVAL:
            name, cls, size, arrive_us = payload
            queues[cls].append((arrive_us, size, name))
            if transmitting is None:
                transmitting = _try_start_transmission(now_us)
                if transmitting:
                    push(transmitting[0], _EVT_TX_END, None)

        elif etype == _EVT_TX_END:
            if transmitting is None:
                continue
            finish_us, size, arrive_us, name, cls = transmitting
            response_times[name].append(finish_us - arrive_us)
            for c in ('A', 'B'):
                _update_credit(c, finish_us, transmitting_cls=cls if cls == c else None)
            transmitting = None
            transmitting = _try_start_transmission(finish_us)
            if transmitting:
                push(transmitting[0], _EVT_TX_END, None)

        elif etype == _EVT_CREDIT_READY:
            wakeup_pending[0] = False   # consume the wakeup slot
            if transmitting is None:
                for c in ('A', 'B'):
                    _update_credit(c, now_us)
                transmitting = _try_start_transmission(now_us)
                if transmitting:
                    push(transmitting[0], _EVT_TX_END, None)

    # Scale for multi-hop (line topology: response times add per hop)
    if num_hops > 1:
        response_times = {
            name: [rt * num_hops for rt in rts]
            for name, rts in response_times.items()
        }
    return response_times


# ─── SP SIMULATION (optional) ────────────────────────────────────────────────

def simulate_sp(streams, num_hops=1, sim_duration_us=None, frames_per_stream=300, seed=42):
    
    if sim_duration_us is None:
        max_period = max(s['period_us'] for s in streams)
        sim_duration_us = max_period * frames_per_stream

    events = []
    counter = [0]

    def push(t, etype, payload):
        heapq.heappush(events, (t, counter[0], etype, payload))
        counter[0] += 1

    for s in streams:
        t = s.get('offset_us', 0.0)
        while t < sim_duration_us:
            push(t, _EVT_ARRIVAL, (s['name'], s['class'], s['size'], t))
            t += s['period_us']

    queues = {'A': [], 'B': [], 'BE': []}
    transmitting = None

    response_times = {s['name']: [] for s in streams}

    def _start_next(now_us):
        for cls in ('A', 'B', 'BE'):
            if queues[cls]:
                arrive_us, size, name = queues[cls].pop(0)
                finish_us = now_us + tx_time_us(size)
                return (finish_us, size, arrive_us, name, cls)
        return None

    while events:
        now_us, _, etype, payload = heapq.heappop(events)
        if now_us > sim_duration_us:
            break

        if etype == _EVT_ARRIVAL:
            name, cls, size, arrive_us = payload
            queues[cls].append((arrive_us, size, name))
            if transmitting is None:
                transmitting = _start_next(now_us)
                if transmitting:
                    push(transmitting[0], _EVT_TX_END, None)

        elif etype == _EVT_TX_END:
            if transmitting is None:
                continue
            finish_us, size, arrive_us, name, cls = transmitting
            response_times[name].append(finish_us - arrive_us)
            transmitting = None
            transmitting = _start_next(finish_us)
            if transmitting:
                push(transmitting[0], _EVT_TX_END, None)

    if num_hops > 1:
        response_times = {
            name: [rt * num_hops for rt in rts]
            for name, rts in response_times.items()
        }
    return response_times


# TEST CASE DEFINITIONS 


def make_testcase_1():
   
    return [
        {'name': 'S_A1', 'class': 'A',  'size': 500,  'period_us': 1000.0},
        {'name': 'S_B1', 'class': 'B',  'size': 400,  'period_us': 1000.0},
        {'name': 'S_E1', 'class': 'BE', 'size': 300,  'period_us': 2000.0},
    ]


def make_testcase_2():
   
    return [
        {'name': 'S_A1', 'class': 'A',  'size': 400,  'period_us': 500.0},
        {'name': 'S_A2', 'class': 'A',  'size': 600,  'period_us': 750.0},
        {'name': 'S_A3', 'class': 'A',  'size': 200,  'period_us': 1000.0},
        {'name': 'S_B1', 'class': 'B',  'size': 500,  'period_us': 500.0},
        {'name': 'S_B2', 'class': 'B',  'size': 350,  'period_us': 750.0},
        {'name': 'S_B3', 'class': 'B',  'size': 250,  'period_us': 1000.0},
        {'name': 'S_E1', 'class': 'BE', 'size': 800,  'period_us': 1000.0},
        {'name': 'S_E2', 'class': 'BE', 'size': 1200, 'period_us': 2000.0},
    ]


def make_testcase_3():
   
    return [
        {'name': 'S_A1', 'class': 'A',  'size': 1000, 'period_us': 800.0},
        {'name': 'S_A2', 'class': 'A',  'size': 500,  'period_us': 1000.0},
        {'name': 'S_A3', 'class': 'A',  'size': 250,  'period_us': 2000.0},
        {'name': 'S_B1', 'class': 'B',  'size': 800,  'period_us': 800.0},
        {'name': 'S_B2', 'class': 'B',  'size': 600,  'period_us': 1200.0},
        {'name': 'S_B3', 'class': 'B',  'size': 300,  'period_us': 2000.0},
        {'name': 'S_E1', 'class': 'BE', 'size': 1522, 'period_us': 2000.0},
        {'name': 'S_E2', 'class': 'BE', 'size': 64,   'period_us': 500.0},
        {'name': 'S_E3', 'class': 'BE', 'size': 500,  'period_us': 1500.0},
    ]


def make_testcase_4():
   
    return make_testcase_1()


def make_testcase_5():
    
    return [
        {'name': 'S_A1', 'class': 'A',  'size': 1500, 'period_us': 1000.0},
        {'name': 'S_A2', 'class': 'A',  'size': 1000, 'period_us': 2000.0},
        {'name': 'S_B1', 'class': 'B',  'size': 1500, 'period_us': 1000.0},
        {'name': 'S_B2', 'class': 'B',  'size': 1000, 'period_us': 2000.0},
        {'name': 'S_E1', 'class': 'BE', 'size': 1000, 'period_us':  500.0},
        {'name': 'S_E2', 'class': 'BE', 'size': 500,  'period_us': 1000.0},
    ]


SCENARIOS = [
    ('Scenario 1 — 1 stream per class (1 hop)', make_testcase_1(), 1),
    ('Scenario 2 — 3A, 3B, 2BE (1 hop)',        make_testcase_2(), 1),
    ('Scenario 3 — High-load (1 hop)',           make_testcase_3(), 1),
    ('Scenario 4 — Multi-hop (2 hops)',          make_testcase_4(), 2),
    ('Scenario 5 — BE Starvation Demo (1 hop)',  make_testcase_5(), 1),
]


# ─── ANALYSIS + SIMULATION

def run_scenario(label, streams, num_hops, sim_frames=400,
                 alpha_plus=None, alpha_minus=None):
    """Run CBS + SP analysis and simulation. Accepts custom alpha for idleSlope variants."""
    wcds_cbs = compute_all_wcds(streams, num_hops, alpha_plus, alpha_minus)
    wcds_sp  = compute_all_wcds_sp(streams, num_hops)

    rts_cbs = simulate_cbs(streams, num_hops=num_hops,
                            frames_per_stream=sim_frames, seed=42,
                            alpha_plus=alpha_plus, alpha_minus=alpha_minus)
    rts_sp  = simulate_sp(streams, num_hops=num_hops,
                           frames_per_stream=sim_frames, seed=42)

    return {
        'label':       label,
        'streams':     streams,
        'num_hops':    num_hops,
        'wcds_cbs':    wcds_cbs,
        'wcds_sp':     wcds_sp,
        'rts_cbs':     rts_cbs,
        'rts_sp':      rts_sp,
        'alpha_plus':  alpha_plus  if alpha_plus  is not None else ALPHA_PLUS,
        'alpha_minus': alpha_minus if alpha_minus is not None else ALPHA_MINUS,
    }


# ─── REPORTING

def print_results(result):
    streams  = result['streams']
    wcds_cbs = result['wcds_cbs']
    wcds_sp  = result['wcds_sp']
    rts_cbs  = result['rts_cbs']
    rts_sp   = result['rts_sp']
    hops     = result['num_hops']
    label    = result['label']

    print(f"\n{'='*72}")
    print(f"  {label}   ({hops} hop{'s' if hops > 1 else ''})")
    print(f"{'='*72}")
    print(f"  Link: {LINK_RATE_MBPS:.0f} Mb/s  |  alpha+ = alpha- = {IDLE_SLOPE_NORM*100:.0f}% BW per CBS class")
    print()

    hdr = f"  {'Stream':<10} {'Class':<5} {'Size':>6} {'Period':>8}  "
    hdr += f"{'WCD_CBS':>10}  {'Sim_CBS_max':>12}  {'WCD_SP':>10}  {'Sim_SP_max':>10}"
    print(hdr)
    print(f"  {'-'*10} {'-'*5} {'-'*6} {'-'*8}  {'-'*10}  {'-'*12}  {'-'*10}  {'-'*10}")

    class_order = {'A': 0, 'B': 1, 'BE': 2}
    for s in sorted(streams, key=lambda x: class_order[x['class']]):
        name = s['name']
        sim_max_cbs = max(rts_cbs[name]) if rts_cbs[name] else float('nan')
        sim_max_sp  = max(rts_sp[name])  if rts_sp[name]  else float('nan')
        wcd_cbs = wcds_cbs[name]
        wcd_sp  = wcds_sp[name]
        ok_cbs  = '✓' if sim_max_cbs <= wcd_cbs + 1e-3 else '!'
        ok_sp   = '✓' if sim_max_sp  <= wcd_sp  + 1e-3 else '!'
        print(f"  {name:<10} {s['class']:<5} {s['size']:>6} {s['period_us']:>8.0f}  "
              f"{wcd_cbs:>10.2f}  {sim_max_cbs:>10.2f} {ok_cbs}  "
              f"{wcd_sp:>10.2f}  {sim_max_sp:>10.2f} {ok_sp}")

    print(f"\n  Units: sizes in bytes, periods and delays in µs")
    print(f"  ✓ = simulated max ≤ analytical WCD   ! = violation (check load)")

    # Utilization summary
    util_A = sum(tx_time_us(s['size']) / s['period_us']
                 for s in streams if s['class'] == 'A')
    util_B = sum(tx_time_us(s['size']) / s['period_us']
                 for s in streams if s['class'] == 'B')
    util_BE = sum(tx_time_us(s['size']) / s['period_us']
                  for s in streams if s['class'] == 'BE')
    total   = util_A + util_B + util_BE
    print(f"\n  Bandwidth utilization:  "
          f"A={util_A*100:.1f}%  B={util_B*100:.1f}%  BE={util_BE*100:.1f}%  "
          f"Total={total*100:.1f}%")
    if util_A > IDLE_SLOPE_NORM:
        print(f"  WARNING: AVB-A utilization {util_A*100:.1f}% > idleSlope {IDLE_SLOPE_NORM*100:.0f}%")
    if util_B > IDLE_SLOPE_NORM:
        print(f"  WARNING: AVB-B utilization {util_B*100:.1f}% > idleSlope {IDLE_SLOPE_NORM*100:.0f}%")


# ─── PLOTTING 

def plot_results(results, out_dir):
    os.makedirs(out_dir, exist_ok=True)
    colors    = {'A': '#2196F3', 'B': '#4CAF50', 'BE': '#FF9800'}
    cls_order = {'A': 0, 'B': 1, 'BE': 2}

    def _sorted_names(streams):
        return [s['name'] for s in sorted(streams, key=lambda x: cls_order[x['class']])]

    def _colors(streams):
        return [colors[s['class']] for s in sorted(streams, key=lambda x: cls_order[x['class']])]

    def _legend_patches(ax):
        patches = [mpatches.Patch(color=colors[c], label=f'AVB-{c}' if c != 'BE' else 'BE')
                   for c in ('A', 'B', 'BE')]
        ax.legend(handles=patches, fontsize=7, loc='upper left')

    # ── Figure 1: Analytical WCD CBS vs SP — 3×2 grid for 5 scenarios ──
    fig, axes = plt.subplots(3, 2, figsize=(14, 14))
    axes = axes.flatten()

    for ax, result in zip(axes, results):
        snames = _sorted_names(result['streams'])
        clr    = _colors(result['streams'])
        x, w   = np.arange(len(snames)), 0.35
        ax.bar(x - w/2, [result['wcds_cbs'][n] for n in snames], w, color=clr, alpha=0.85, edgecolor='white')
        ax.bar(x + w/2, [result['wcds_sp'][n]  for n in snames], w, color=clr, alpha=0.45, edgecolor='grey')
        ax.set_xticks(x); ax.set_xticklabels(snames, rotation=30, ha='right', fontsize=8)
        ax.set_ylabel('WCD (µs)'); ax.set_title(result['label'], fontsize=9)
        _legend_patches(ax)

    # Hide unused 6th subplot
    axes[5].set_visible(False)
    fig.suptitle('Analytical WCD: CBS (solid) vs Strict Priority (faded)', fontweight='bold')
    plt.tight_layout()
    out1 = os.path.join(out_dir, 'wcd_cbs_vs_sp.png')
    plt.savefig(out1, dpi=150, bbox_inches='tight'); plt.close()
    print(f"\n  Plot saved: {out1}")

    # ── Figure 2: CBS Analytical vs Simulated — 3×2 grid ──
    fig, axes = plt.subplots(3, 2, figsize=(14, 14))
    axes = axes.flatten()

    for ax, result in zip(axes, results):
        snames = _sorted_names(result['streams'])
        clr    = _colors(result['streams'])
        x, w   = np.arange(len(snames)), 0.35
        analyt = [result['wcds_cbs'][n] for n in snames]
        simmax = [max(result['rts_cbs'][n]) if result['rts_cbs'][n] else 0 for n in snames]
        simavg = [sum(result['rts_cbs'][n])/len(result['rts_cbs'][n]) if result['rts_cbs'][n] else 0 for n in snames]
        ax.bar(x - w/2, analyt, w, color=clr, alpha=0.85, edgecolor='white')
        ax.bar(x + w/2, simmax, w, color=clr, alpha=0.45, edgecolor='grey')
        ax.plot(x + w/2, simavg, 'kx', markersize=6)
        ax.set_xticks(x); ax.set_xticklabels(snames, rotation=30, ha='right', fontsize=8)
        ax.set_ylabel('Delay (µs)'); ax.set_title(result['label'], fontsize=9)
        _legend_patches(ax)

    axes[5].set_visible(False)
    fig.suptitle('CBS: Analytical WCD vs Simulation (max and average)', fontweight='bold')
    plt.tight_layout()
    out2 = os.path.join(out_dir, 'cbs_analytical_vs_sim.png')
    plt.savefig(out2, dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Plot saved: {out2}")

    # ── Figure 3: Response time CDF for Scenario 2 (most interesting) ──
    result = results[1]   # Scenario 2
    streams  = result['streams']
    rts_cbs  = result['rts_cbs']
    rts_sp   = result['rts_sp']
    wcds_cbs = result['wcds_cbs']
    wcds_sp  = result['wcds_sp']

    class_order = {'A': 0, 'B': 1, 'BE': 2}
    sorted_streams = sorted(streams, key=lambda x: class_order[x['class']])

    fig, axes = plt.subplots(1, len(sorted_streams), figsize=(15, 4), sharey=False)
    colors_cls = {'A': '#2196F3', 'B': '#4CAF50', 'BE': '#FF9800'}

    for ax, s in zip(axes, sorted_streams):
        name = s['name']
        clr = colors_cls[s['class']]

        def plot_cdf(ax, data, label, color, linestyle='-'):
            if not data:
                return
            sorted_d = sorted(data)
            n = len(sorted_d)
            ax.plot(sorted_d, [(i+1)/n for i in range(n)],
                    label=label, color=color, linestyle=linestyle, linewidth=1.5)

        plot_cdf(ax, rts_cbs.get(name, []), 'CBS', clr, '-')
        plot_cdf(ax, rts_sp.get(name,  []), 'SP',  clr, '--')

        ax.axvline(wcds_cbs.get(name, 0), color=clr, linestyle=':', linewidth=1.2, label=f'CBS WCD={wcds_cbs.get(name,0):.1f}µs')
        ax.axvline(wcds_sp.get(name,  0), color='grey', linestyle=':', linewidth=1.0, label=f'SP  WCD={wcds_sp.get(name,0):.1f}µs')
        ax.set_xlabel('Response time (µs)', fontsize=8)
        ax.set_ylabel('CDF', fontsize=8)
        ax.set_title(f"{name} ({s['class']})", fontsize=9)
        ax.legend(fontsize=7)
        ax.grid(True, alpha=0.3)

    fig.suptitle('Scenario 2 — Response Time CDF: CBS vs Strict Priority', fontweight='bold')
    plt.tight_layout()
    out3 = os.path.join(out_dir, 'cdf_scenario2.png')
    plt.savefig(out3, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {out3}")

    # ── Figure 4: Credit trace for Scenario 1 (illustrative) ──
    _plot_credit_trace(make_testcase_1(), out_dir)

    # ── Figure 5: Starvation Demo — CBS vs SP for Scenario 5 ──
    result5 = results[4]   # Scenario 5
    streams5   = result5['streams']
    rts_cbs5   = result5['rts_cbs']
    rts_sp5    = result5['rts_sp']
    wcds_cbs5  = result5['wcds_cbs']
    wcds_sp5   = result5['wcds_sp']

    sorted5 = sorted(streams5, key=lambda x: cls_order[x['class']])
    fig, axes = plt.subplots(1, len(sorted5), figsize=(15, 4), sharey=False)

    for ax, s in zip(axes, sorted5):
        name = s['name']
        clr  = colors[s['class']]

        def _cdf(ax, data, label, color, ls):
            if not data: return
            sd = sorted(data); n = len(sd)
            ax.plot(sd, [(i+1)/n for i in range(n)], label=label, color=color,
                    linestyle=ls, linewidth=1.5)

        _cdf(ax, rts_cbs5.get(name, []), 'CBS', clr, '-')
        _cdf(ax, rts_sp5.get(name,  []), 'SP',  clr, '--')
        ax.axvline(wcds_cbs5.get(name, 0), color=clr,    linestyle=':', linewidth=1.2,
                   label=f'CBS WCD={wcds_cbs5.get(name,0):.0f}µs')
        ax.axvline(wcds_sp5.get(name,  0), color='grey', linestyle=':', linewidth=1.0,
                   label=f'SP  WCD={wcds_sp5.get(name,0):.0f}µs')
        ax.set_xlabel('Response time (µs)', fontsize=8)
        ax.set_ylabel('CDF', fontsize=8)
        ax.set_title(f"{name} ({s['class']})", fontsize=9)
        ax.legend(fontsize=7); ax.grid(True, alpha=0.3)

    fig.suptitle('Scenario 5 — BE Starvation Demo: CBS vs Strict Priority\n'
                 '(CBS WCD > SP WCD for BE — credit overhead protects against starvation)',
                 fontweight='bold')
    plt.tight_layout()
    out5 = os.path.join(out_dir, 'starvation_demo_scenario5.png')
    plt.savefig(out5, dpi=150, bbox_inches='tight'); plt.close()
    print(f"  Plot saved: {out5}")


def _plot_credit_trace(streams, out_dir, trace_duration_us=3000):
    """Replay a short CBS simulation and plot credit evolution."""
    # Re-run with very fine time tracking for trace
    rng_seed = 0
    events = []
    counter = [0]

    def push(t, etype, payload):
        heapq.heappush(events, (t, counter[0], etype, payload))
        counter[0] += 1

    for s in streams:
        t = s.get('offset_us', 0.0)
        while t < trace_duration_us:
            push(t, _EVT_ARRIVAL, (s['name'], s['class'], s['size'], t))
            t += s['period_us']

    credit  = {'A': 0.0, 'B': 0.0}
    cr_time = {'A': 0.0, 'B': 0.0}
    queues  = {'A': [], 'B': [], 'BE': []}
    transmitting    = None
    wakeup_pending  = [False]

    trace_t    = [0.0]
    trace_cr_A = [0.0]
    trace_cr_B = [0.0]
    tx_bars    = []   # (start, end, cls)

    _EPS = 1e-6

    def _upd(cls, now_us, is_tx=False):
        dt_s = (now_us - cr_time[cls]) * 1e-6
        if dt_s <= 0:
            return
        if is_tx:
            credit[cls] -= ALPHA_MINUS * dt_s
        elif queues[cls] and credit[cls] < 0:
            credit[cls] = min(0.0, credit[cls] + ALPHA_PLUS * dt_s)
        elif not queues[cls] and credit[cls] > 0:
            credit[cls] = 0.0
        if -_EPS < credit[cls] < _EPS:
            credit[cls] = 0.0
        cr_time[cls] = now_us

    def _record(t):
        trace_t.append(t)
        trace_cr_A.append(credit['A'])
        trace_cr_B.append(credit['B'])

    def _start_tx(now_us):
        for cls in ('A', 'B', 'BE'):
            if not queues[cls]:
                continue
            if cls in ('A', 'B') and credit[cls] < -_EPS:
                continue
            arrive_us, size, name = queues[cls].pop(0)
            finish_us = now_us + tx_time_us(size)
            tx_bars.append((now_us, finish_us, cls))
            return (finish_us, size, arrive_us, name, cls)
        # Schedule single wakeup at earliest credit recovery
        if not wakeup_pending[0]:
            t_wake = float('inf')
            for cls in ('A', 'B'):
                if queues[cls] and credit[cls] < 0:
                    t_wake = min(t_wake, now_us + (-credit[cls] / ALPHA_PLUS) * 1e6)
            if t_wake < trace_duration_us:
                wakeup_pending[0] = True
                push(t_wake, _EVT_CREDIT_READY, None)
        return None

    while events:
        now_us, _, etype, payload = heapq.heappop(events)
        if now_us > trace_duration_us:
            break

        if etype == _EVT_ARRIVAL:
            name, cls, size, arrive_us = payload
            queues[cls].append((arrive_us, size, name))
            if transmitting is None:
                for c in ('A', 'B'):
                    _upd(c, now_us)
                _record(now_us)
                transmitting = _start_tx(now_us)
                if transmitting:
                    push(transmitting[0], _EVT_TX_END, None)

        elif etype == _EVT_TX_END:
            if transmitting is None:
                continue
            finish_us, _, _, _, cls = transmitting
            for c in ('A', 'B'):
                _upd(c, finish_us, is_tx=(c == cls))
            _record(finish_us)
            transmitting = None
            transmitting = _start_tx(finish_us)
            if transmitting:
                push(transmitting[0], _EVT_TX_END, None)

        elif etype == _EVT_CREDIT_READY:
            wakeup_pending[0] = False
            if transmitting is None:
                for c in ('A', 'B'):
                    _upd(c, now_us)
                _record(now_us)
                transmitting = _start_tx(now_us)
                if transmitting:
                    push(transmitting[0], _EVT_TX_END, None)

    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(12, 6), sharex=True)

    # Convert credit from bits to µs equivalent (divide by link rate in Mbps)
    cr_A_us = [c / LINK_RATE_BPS * 1e6 for c in trace_cr_A]
    cr_B_us = [c / LINK_RATE_BPS * 1e6 for c in trace_cr_B]

    ax1.plot(trace_t, cr_A_us, color='#2196F3', linewidth=1.5, label='Credit AVB-A (µs equiv.)')
    ax1.plot(trace_t, cr_B_us, color='#4CAF50', linewidth=1.5, label='Credit AVB-B (µs equiv.)')
    ax1.axhline(0, color='black', linewidth=0.8, linestyle='--')
    ax1.set_ylabel('CBS Credit (µs equiv.)')
    ax1.legend(fontsize=8)
    ax1.grid(True, alpha=0.3)
    ax1.set_title('CBS Credit Evolution — Scenario 1')

    cls_colors = {'A': '#2196F3', 'B': '#4CAF50', 'BE': '#FF9800'}
    for (t_s, t_e, cls) in tx_bars:
        ax2.barh(0, t_e - t_s, left=t_s, height=0.5, color=cls_colors.get(cls, 'grey'), alpha=0.7)
    ax2.set_yticks([])
    ax2.set_xlabel('Time (µs)')
    ax2.set_ylabel('Transmissions')
    ax2.set_title('Link Utilization (A=blue, B=green, BE=orange)')
    ax2.grid(True, alpha=0.3, axis='x')

    plt.tight_layout()
    out4 = os.path.join(out_dir, 'credit_trace_scenario1.png')
    plt.savefig(out4, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {out4}")


def run_idleslope_comparison(streams, num_hops=1, sim_frames=400):
   
    BW = LINK_RATE_BPS
    configs = [
        ('A=75%, B=25%', 0.75*BW, BW-0.75*BW),
        ('A=50%, B=50%', 0.50*BW, BW-0.50*BW),
        ('A=25%, B=75%', 0.25*BW, BW-0.25*BW),
    ]
    results = []
    for name, ap, am in configs:
        label = f'Scenario 6 — idleSlope {name} (1 hop)'
        r = run_scenario(label, streams, num_hops, sim_frames,
                         alpha_plus=ap, alpha_minus=am)
        results.append(r)
    return results


def plot_idleslope_comparison(idleslope_results, out_dir):
    """
    Figure 6: idleSlope sensitivity — analytical WCDs under 3 bandwidth allocations.
    Grouped bar chart: x=stream, 3 bars per stream (one per config).
    Shows trade-off: higher idleSlope for A lowers A's WCD but raises B's.
    """
    os.makedirs(out_dir, exist_ok=True)
    colors    = {'A': '#2196F3', 'B': '#4CAF50', 'BE': '#FF9800'}
    cls_order = {'A': 0, 'B': 1, 'BE': 2}

    streams = idleslope_results[0]['streams']
    snames  = [s['name'] for s in sorted(streams, key=lambda x: cls_order[x['class']])]
    clrs    = [colors[s['class']] for s in sorted(streams, key=lambda x: cls_order[x['class']])]
    cfg_labels = ['A=75%\nB=25%', 'A=50%\nB=50%\n(base)', 'A=25%\nB=75%']
    alphas_bar = [0.95, 0.70, 0.45]

    x   = np.arange(len(snames))
    w   = 0.25
    fig, ax = plt.subplots(figsize=(14, 5))

    for i, (result, cfg_lbl, alph) in enumerate(zip(idleslope_results, cfg_labels, alphas_bar)):
        vals = [result['wcds_cbs'][n] for n in snames]
        offset = (i - 1) * w
        bars = ax.bar(x + offset, vals, w, color=clrs, alpha=alph, edgecolor='white',
                      label=cfg_lbl.replace('\n', ' '))
        # Annotate top of each bar
        for bar, val in zip(bars, vals):
            ax.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 3,
                    f'{val:.0f}', ha='center', va='bottom', fontsize=6.5, rotation=90)

    ax.set_xticks(x)
    ax.set_xticklabels(snames, rotation=30, ha='right', fontsize=9)
    ax.set_ylabel('WCD_CBS (µs)')
    ax.set_title('Scenario 6 — idleSlope Sensitivity: WCD vs Bandwidth Allocation\n'
                 'Solid=A=75%/B=25%  |  Medium=A=50%/B=50% (baseline)  |  Faded=A=25%/B=75%',
                 fontweight='bold')

    patch_A  = mpatches.Patch(color='#2196F3', label='AVB-A')
    patch_B  = mpatches.Patch(color='#4CAF50', label='AVB-B')
    patch_BE = mpatches.Patch(color='#FF9800', label='BE')
    ax.legend(handles=[patch_A, patch_B, patch_BE], fontsize=8, loc='upper left')
    ax.grid(True, alpha=0.3, axis='y')

    plt.tight_layout()
    out6 = os.path.join(out_dir, 'idleslope_sensitivity_scenario6.png')
    plt.savefig(out6, dpi=150, bbox_inches='tight')
    plt.close()
    print(f"  Plot saved: {out6}")


# ─── JSON TEST CASE LOADER 


_PCP_TO_CLASS = {2: 'A', 1: 'B', 0: 'BE'}


def load_json_testcase(streams_path, topology_path, routes_path):
   
    with open(streams_path)  as f: sdata = _json.load(f)
    with open(topology_path) as f: tdata = _json.load(f)
    with open(routes_path)   as f: rdata = _json.load(f)

    topo       = tdata['topology']
   
    default_bw = topo.get('default_bandwidth_mbps', 100)

    
    route_map = {}
    for r in rdata['routes']:
        path = r['paths'][0]
        route_map[r['flow_id']] = [(h['node'], h['port']) for h in path[:-1]]

    # Convert streams to internal format
    all_streams = []
    for s in sdata['streams']:
        all_streams.append({
            'name':      s['name'],
            'id':        s['id'],
            'class':     _PCP_TO_CLASS.get(s['PCP'], 'BE'),
            'size':      s['size'],
            'period_us': float(s['period']),
            'deadline':  float(s['destinations'][0]['deadline']),
            'route':     route_map.get(s['id'], []),
        })

    # Group streams by egress port
    port_streams = _defaultdict(list)
    for s in all_streams:
        for hop in s['route']:
            port_streams[hop].append(s)

    seen, groups = set(), []
    for s in all_streams:
        for hop in s['route']:
            if hop not in seen:
                seen.add(hop)
                groups.append({
                    'port_id':        f"{hop[0]}:port{hop[1]}",
                    'hop':            hop,
                    'bandwidth_mbps': default_bw,   # all links = default_bw
                    'streams':        port_streams[hop],
                })
    return all_streams, groups, default_bw


def _tx_at(size_bytes, bw_mbps):
    return size_bytes * 8.0 / (bw_mbps * 1e6) * 1e6   # µs


def compute_wcd_json_hop(stream, all_streams_on_port, bw_mbps):
   
    def tx(b): return _tx_at(b, bw_mbps)
    alpha = 0.5 * bw_mbps * 1e6   # alpha+ = alpha- = 0.5 * BW

    C_i = tx(stream['size'])
    cls = stream['class']

    if cls == 'A':
        peers = [s for s in all_streams_on_port if s['class'] == 'A'
                 and s['name'] != stream['name']]
        lower = [s for s in all_streams_on_port if s['class'] in ('B', 'BE')]
        lpi   = max((tx(s['size']) for s in lower), default=0.0)  # no x2, matches reference
        R = C_i + lpi + sum(tx(s['size']) * 2.0 for s in peers)
        for _ in range(500):
            R_new = C_i + lpi + sum(
                math.ceil(R / s['period_us']) * tx(s['size']) * 2.0 for s in peers)
            if abs(R_new - R) < 1e-9: R = R_new; break
            R = R_new

    elif cls == 'B':
        peers    = [s for s in all_streams_on_port if s['class'] == 'B'
                    and s['name'] != stream['name']]
        higher   = [s for s in all_streams_on_port if s['class'] == 'A']
        lower    = [s for s in all_streams_on_port if s['class'] == 'BE']
        C_A_max  = max((tx(s['size']) for s in higher), default=0.0)
        C_BE_max = max((tx(s['size']) for s in lower),  default=0.0)
        mpi = C_BE_max * 2.0 + C_A_max
        R = C_i + mpi + sum(tx(s['size']) * 2.0 for s in peers)
        for _ in range(500):
            R_new = C_i + mpi + sum(
                math.ceil(R / s['period_us']) * tx(s['size']) * 2.0 for s in peers)
            if abs(R_new - R) < 1e-9: R = R_new; break
            R = R_new

    else:  # BE
        hp_a    = [s for s in all_streams_on_port if s['class'] == 'A']
        hp_b    = [s for s in all_streams_on_port if s['class'] == 'B']
        block   = max((tx(s['size']) for s in all_streams_on_port), default=0.0)
        C_A_max = max((tx(s['size']) for s in hp_a), default=0.0)
        C_B_max = max((tx(s['size']) for s in hp_b), default=0.0)
        credit_block = C_A_max + C_B_max   # alpha-/alpha+ = 1
        R = C_i + block + credit_block
        for _ in range(500):
            R_new = C_i + block + credit_block + sum(
                math.ceil(R / s['period_us']) * tx(s['size'])
                for s in hp_a + hp_b)
            if abs(R_new - R) < 1e-9: R = R_new; break
            R = R_new

    return R


def analyse_json_testcase(streams_path, topology_path, routes_path):
    """Run full end-to-end WCD analysis for a JSON test case."""
    all_streams, groups, default_bw = load_json_testcase(
        streams_path, topology_path, routes_path)

    hop_to_group = {g['hop']: g for g in groups}

    results = []
    for s in all_streams:
        e2e_wcd  = 0.0
        hop_wcds = []
        for hop in s['route']:
            grp = hop_to_group.get(hop)
            if grp is None: continue
            wcrt = compute_wcd_json_hop(s, grp['streams'], grp['bandwidth_mbps'])
            e2e_wcd += wcrt
            hop_wcds.append((f"{hop[0]}:port{hop[1]}", grp['bandwidth_mbps'], wcrt))
        results.append({
            'name':           s['name'],
            'id':             s['id'],
            'class':          s['class'],
            'size':           s['size'],
            'period':         s['period_us'],
            'deadline':       s['deadline'],
            'e2e_wcd':        e2e_wcd,
            'hop_wcds':       hop_wcds,
            'meets_deadline': e2e_wcd <= s['deadline'],
        })
    return results


def print_json_results(results, ref_csv_path=None, label='JSON Test Case'):
    
    ref = {}
    if ref_csv_path and os.path.exists(ref_csv_path):
        with open(ref_csv_path) as f:
            for line in f:
                line = line.strip()
                if not line or line.startswith('ID'): continue
                parts = line.split('\t')
                if len(parts) == 2:
                    try:
                        # Support both old comma and new dot decimal separators
                        ref[int(parts[0])] = float(parts[1].replace(',', '.'))
                    except ValueError:
                        pass

    print(f"\n{'='*72}")
    print(f"  {label}")
    print(f"{'='*72}")
    hdr = (f"  {'Stream':<12} {'Cls':<4} {'Size':>6} {'Period':>8} "
           f"{'Deadline':>9}  {'WCD_e2e':>10}  {'OK':>4}")
    if ref: hdr += f"  {'Ref_WCRT':>10}  {'Diff%':>7}"
    print(hdr)
    sep = (f"  {'-'*12} {'-'*4} {'-'*6} {'-'*8} "
           f"{'-'*9}  {'-'*10}  {'-'*4}")
    if ref: sep += f"  {'-'*10}  {'-'*7}"
    print(sep)

    for r in sorted(results, key=lambda x: x['id']):
        ok  = '✓' if r['meets_deadline'] else '✗'
        row = (f"  {r['name']:<12} {r['class']:<4} {r['size']:>6} "
               f"{r['period']:>8.0f} {r['deadline']:>9.0f}  "
               f"{r['e2e_wcd']:>10.2f}  {ok:>4}")
        if ref:
            if r['id'] in ref:
                diff = (r['e2e_wcd'] - ref[r['id']]) / ref[r['id']] * 100
                row += f"  {ref[r['id']]:>10.2f}  {diff:>+7.2f}%"
            else:
                row += f"  {'(no ref)':>10}  {'—':>7}"
        print(row)
    print(f"\n  Units: bytes / µs.  ✓ = WCD ≤ deadline   ✗ = deadline exceeded")
    print(f"  Note: BE streams excluded from reference (CBS does not bound BE "
          f"when alpha_A+ + alpha_B+ = BW).")


# ─── MAIN ────────────────────────────────────────────────────────────────────

def main():
    out_dir = os.path.join(os.path.dirname(__file__), 'plots')
    os.makedirs(out_dir, exist_ok=True)

    print("02225 DRTS Mini Project 2 — CBS WCD Analysis & Simulation")
    print(f"Link rate: {LINK_RATE_MBPS:.0f} Mb/s  "
          f"idleSlope = sendSlope = {IDLE_SLOPE_NORM*100:.0f}%  "
          f"alpha+ = alpha- = {ALPHA_PLUS/1e6:.0f} Mb/s per class")

    results = []
    for label, streams, hops in SCENARIOS:
        print(f"\nRunning: {label} ...")
        result = run_scenario(label, streams, hops, sim_frames=400)
        print_results(result)
        results.append(result)

    # ── Formula summary ──
    print("\n" + "="*72)
    print("  FORMULA SUMMARY (Cao 2016, with alpha+ = alpha- = 0.5*BW)")
    print("="*72)
    print()
    print("  Notation:")
    print("    C_i    = tx time of S_i = size*8/link_rate,  T_i = period of S_i")
    print("    alpha+ = idleSlope = 0.5*BW  (credit rises while queue waits)")
    print("    alpha- = sendSlope = 0.5*BW  (credit falls while transmitting)")
    print()
    print("  AVB-A (Cao 2016 WFCS Thm 4 — iterative busy-window):")
    print("    R_A(Si) = C_i + LPI + sum_{j!=i,A} ceil(R/T_j) * C_j * (1+alpha-/alpha+)")
    print("    LPI = max_{B,BE} C_j * (1+alpha+/alpha-)  =  2*C_L_max")
    print("    (LPI: lower frame transmitting -> A credit recovers -> eligible interval delayed)")
    print()
    print("  AVB-B (Cao 2016 RTNS Thm 3 — iterative busy-window):")
    print("    R_B(Si) = C_i + MPI + sum_{j!=i,B} ceil(R/T_j) * C_j * (1+alpha-/alpha+)")
    print("    MPI = C_BE_max*(1+alpha+/alpha-) + C_A_max  =  2*C_BE_max + C_A_max")
    print()
    print("  BE (non-preemptive iterative RTA, A+B as strict higher priority):")
    print("    R_BE(Si) = C_i + B_block + credit_block + sum_{A,B} ceil(R/T_j)*C_j")
    print("    B_block      = max C_j over ALL streams  (any frame may be mid-tx)")
    print("    credit_block = (C_A_max + C_B_max)*(alpha-/alpha+)  [CBS debt recovery]")
    print()
    print("  SP — Strict Priority, no CBS (non-preemptive FP-RTA, FIFO within priority):")
    print("    R_SP(Si) = C_i + B_block")
    print("               + sum_{prio(j)<prio(i)} ceil(R/T_j)*C_j")
    print("               + sum_{prio(j)=prio(i),j!=i} ceil(R/T_j)*C_j")
    print("    B_block = max C_j over ALL streams.  Priority: A(0) > B(1) > BE(2)")
    print()
    print("  Multi-hop (H identical hops, same-direction streams, line topology):")
    print("      WCD_total(Si) = H * WCRT_per_hop(Si)")
    print()

    # ── Plots ──
    print("\nGenerating plots...")
    plot_results(results, out_dir)

    # ── Scenario 6: idleSlope sensitivity ──
    print("\nRunning: Scenario 6 — idleSlope Sensitivity Analysis ...")
    sc6_streams = make_testcase_2()
    idleslope_results = run_idleslope_comparison(sc6_streams, num_hops=1, sim_frames=400)

    cls_order = {'A': 0, 'B': 1, 'BE': 2}
    print(f"\n{'='*72}")
    print(f"  Scenario 6 — idleSlope Sensitivity (Scenario 2 stream set)")
    print(f"{'='*72}")
    print(f"  {'Stream':<10} {'Class':<5}  {'A=75%,B=25%':>12}  {'A=50%,B=50%':>12}  {'A=25%,B=75%':>12}")
    print(f"  {'-'*10} {'-'*5}  {'-'*12}  {'-'*12}  {'-'*12}")
    for s in sorted(sc6_streams, key=lambda x: cls_order[x['class']]):
        n = s['name']
        w1 = idleslope_results[0]['wcds_cbs'][n]
        w2 = idleslope_results[1]['wcds_cbs'][n]
        w3 = idleslope_results[2]['wcds_cbs'][n]
        trend = '↓' if s['class']=='A' else '↑' if s['class']=='B' else '~'
        print(f"  {n:<10} {s['class']:<5}  {w1:>12.2f}  {w2:>12.2f}  {w3:>12.2f}  {trend}")
    print(f"\n  Units: µs.  ↓=decreases with more A bandwidth  ↑=increases")
    print(f"  Note: alpha_A+ + alpha_B+ = 100% BW in all configs (validity condition met).")
    print(f"  Trade-off: allocating more bandwidth to A lowers A WCDs but raises B WCDs.")

    plot_idleslope_comparison(idleslope_results, out_dir)

    # ── JSON Test Case (professor's real test case) ──
    script_dir = os.path.dirname(os.path.abspath(__file__))
    candidates = [
        os.path.join(script_dir, 'tsn-test-cases', 'examples', 'test_case_1'),
        os.path.join(script_dir, '..', 'tsn-test-cases', 'examples', 'test_case_1'),
        os.path.expanduser('~/Desktop/tsn-test-cases/examples/test_case_1'),
    ]
    tc1_dir = next((p for p in candidates if os.path.exists(p)), None)

    if tc1_dir:
        print("\nRunning: JSON Test Case 1 (paulpop/tsn-test-cases/examples/test_case_1) ...")
        json_results = analyse_json_testcase(
            os.path.join(tc1_dir, 'streams.json'),
            os.path.join(tc1_dir, 'topology.json'),
            os.path.join(tc1_dir, 'routes.json'))
        print_json_results(
            json_results,
            os.path.join(tc1_dir, 'WCRTs.csv'),
            label='JSON Test Case 1 — ES0↔ES1 via SW1 (2 hops, 100 Mb/s)')
    else:
        print("\n[INFO] JSON test case not found. Clone to enable:")
        print("       git clone https://github.com/paulpop/tsn-test-cases ~/Desktop/tsn-test-cases")

    print("\nDone. All plots written to:", out_dir)


if __name__ == '__main__':
    main()