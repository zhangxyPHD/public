import pandas as pd
import matplotlib.pyplot as plt
import os
import numpy as np
import matplotlib.ticker as ticker

# Wes=("1" "5" "10")
# Ohs=("0.0035" "0.01")
# J1s="0 0.01 0.05 0.1 0.15"
# J2s="0.2 0.25 0.3 0.35 0.4"

# Configuration parameters
# variables = [0,0.005,0.01,0.025,0.05,0.075,0.1,0.125,0.15,0.175,0.2,0.4]
variables = [0,0.01,0.05,0.1,0.15,0.2,0.25,0.3,0.35,0.4,0.6,0.8,1.0,2.0,4.0,6.0,8.0,10.0]
Oh=0.01
We=5
variables_exist=[]
variable_name="J"
currentDir=os.path.dirname(os.path.abspath(__file__))
results_folder = os.path.join(currentDir, "Results_Running")
plot_folder = os.path.join(currentDir, f"We{We}_J_Oh{Oh}")
if not os.path.exists(plot_folder):
    os.makedirs(plot_folder)
# Data storage dictionary
data = {}
plt.rc('font', size=14)          # Default text size
plt.rc('axes', titlesize=18)     # Title font size
plt.rc('axes', labelsize=16)     # X and Y label size
plt.rc('xtick', labelsize=14)    # X tick label size
plt.rc('ytick', labelsize=14)    # Y tick label size
plt.rc('legend', fontsize=14)    # Legend font size
plt.rc('figure', titlesize=20)   # Figure title font size
# Load data
for variable in variables:
    J=variable
    maxLevel=11
    epsilon=0.001
    filename = f"Bo0.0-We{We}-J{J}-Oh{Oh}-MAXlevel{maxLevel}-epsilon{epsilon}.csv"
    filename1 = f"Bo0.0-We{We}-J{J}-Oh{Oh}-MAXlevel10-epsilon{epsilon}.csv"
    filepath = os.path.join(results_folder, filename)
    filepath1 = os.path.join(results_folder, filename1)
    if not os.path.exists(filepath):
        print(f"{filename} does not exist")        
        if not os.path.exists(filepath1):
            print(f"{filename1} does not exist")
            continue
        else:
            df = pd.read_csv(filepath1)
    else:
        df = pd.read_csv(filepath)
    variables_exist.append(variable)
    t = df["t"]
    vcm1 = df["vcm1"]
    dvcm1_dt = np.gradient(vcm1, t)
    dvcm1_dt[0] = 0
    uH = df["uH"]
    vR = df["vR"]
    ke = df["ke"]

    mask1 = (t >= 0) & (t <= 1)
    dvcm1_dt_sub1 = dvcm1_dt[mask1]
    t_sub1 = t[mask1]
    if len(dvcm1_dt_sub1) > 0:
        F1 = np.max(dvcm1_dt_sub1)
        F1_time = t_sub1[np.argmax(dvcm1_dt_sub1)]
        mask2 = (t >= 1.5) & (t <= 5)

    dvcm1_dt_sub2 = dvcm1_dt[mask2]
    t_sub2 = t[mask2].values
    if len(dvcm1_dt_sub2) > 0:
        F2 = np.max(dvcm1_dt_sub2)
        F2_time = t_sub2[np.argmax(dvcm1_dt_sub2)]

    Zmax = df["Zmax"].values
    if len(uH) > 0:
        uHMAX = np.max(uH)
        uHMAX_t = t[np.argmax(uH)]

    Rmax = df["Rmax"].values
    if len(Rmax) > 0:
        RmaxMAX = np.max(Rmax)
        RmaxMAX_t = t[np.argmax(Rmax)]

    if len(vR) > 0:
        vRMAX = np.max(vR)
        vRMAX_t = t[np.argmax(vR)]

    data[variable] = {"t": t, "vcm1": vcm1, "dvcm1_dt": dvcm1_dt, "uH": uH, "Zmax":Zmax, "vR": vR, "Rmax":Rmax, "F1_t":[F1,F1_time], "F2_t":[F2,F2_time], "uHMAX_t":[uHMAX,uHMAX_t], "vRMAX_t":[vRMAX,vRMAX_t], "RmaxMAX_t":[RmaxMAX,RmaxMAX_t], "ke":ke}

# Plotting function
def plot_figure(x, y, labels, xlabel, ylabel, title, save_as, xlim=None, ylim=None, xticks=None, yticks=None,xlog=False, ylog=False):
    plt.figure(figsize=(10,6))
    for label, (x_vals, y_vals) in labels.items():
        plt.plot(x_vals, y_vals, label=f"{variable_name} = {label}")
    plt.legend()
    plt.xlabel(xlabel)
    plt.ylabel(ylabel)
    plt.title(title)
    if xlim:
        plt.xlim(xlim)
    if ylim:
        plt.ylim(ylim)
    if xticks:
        plt.xticks(np.arange(*xticks))
    if yticks:
        plt.yticks(np.arange(*yticks))
    if xlog:
        plt.xscale('log')
    if ylog:
        plt.yscale('log')

    plt.legend(loc='upper left', bbox_to_anchor=(1, 1))
    # 调整布局以防止图例遮挡图形
    plt.tight_layout(rect=[0, 0, 0.85, 1])
    plt.savefig(save_as, dpi=300, bbox_inches='tight')
    plt.close()

# Plot mv-t
labels_mv = {variable:(data[variable]["t"], data[variable]["vcm1"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_mv,
    xlabel="t", ylabel="mv",
    title=f"mv-t for Different {variable_name}",
    save_as=os.path.join(plot_folder, f"figure_mv_vs_t_for_different_{variable_name}.png")
)

# Plot F-t
labels_F = {variable:(data[variable]["t"], data[variable]["dvcm1_dt"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_F,
    xlabel="t", ylabel="d(mv)/dt",
    title=f"F-t for Different {variable_name}",
    # xlim=(-0.5, 5), ylim=(0, 4.5),
    save_as=os.path.join(plot_folder, f"figure_F_vs_t_for_different_{variable_name}.png")
)

# Plot inset F-t
plot_figure(
    x=None, y=None, labels=labels_F,
    xlabel="t", ylabel="d(mv)/dt",
    title=f"Inset F-t for Different {variable_name}",
    save_as=os.path.join(plot_folder, f"figure_inset_F_vs_t_for_different_{variable_name}.png"),
    xlim=(3.0, 3.8), ylim=(-1, 5),
    xticks=(3.0, 3.8, 0.2), yticks=(-1, 5, 1)
)

# Plot ZMAXimax-t
labels_ZMAX = {variable:(data[variable]["t"], data[variable]["Zmax"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_ZMAX,
    xlabel="t", ylabel="ZmaxAxi",
    title=f"ZmaxAxi-t for Different {variable_name}",
    ylog=True,
    save_as=os.path.join(plot_folder, f"figure_ZmaxAxi_vs_t_for_different_{variable_name}.png")
)

# Plot uH-t
labels_uH = {variable:(data[variable]["t"], data[variable]["uH"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_uH,
    xlabel="t", ylabel="uH",
    title=f"uH-t for Different {variable_name}",
    # ylim=(-2, 20),
    save_as=os.path.join(plot_folder, f"figure_uH_vs_t_for_different_{variable_name}.png")
)

# Plot vR-t
labels_vR = {variable:(data[variable]["t"], data[variable]["vR"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_vR,
    xlabel="t", ylabel="vR",
    title=f"vR-t for Different {variable_name}",
    # ylim=(-2, 20),
    save_as=os.path.join(plot_folder, f"figure_vR_vs_t_for_different_{variable_name}.png")
)
# Plot ke-t
labels_ke = {variable:(data[variable]["t"], data[variable]["ke"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_ke,
    xlabel="t", ylabel="ke",
    title=f"ke-t for Different {variable_name}",
    save_as=os.path.join(plot_folder, f"figure_ke_vs_t_for_different_{variable_name}.png")
)
# Plot Rmax-t
labels_Rmax = {variable:(data[variable]["t"], data[variable]["Rmax"]) for variable in data}
plot_figure(
    x=None, y=None, labels=labels_Rmax,
    xlabel="t", ylabel="Rmax",
    title=f"Rmax-t for Different {variable_name}",
    # ylim=(-2, 20),
    save_as=os.path.join(plot_folder, f"figure_Rmax_vs_t_for_different_{variable_name}.png")
)
####################################################################################################
# Plot F1, F2 vs variable
F1_values = [data[variable]["F1_t"][0] for variable in data]
F2_values = [data[variable]["F2_t"][0] for variable in data]
plt.figure(figsize=(8, 6))
plt.plot(variables_exist, F1_values, 'o-', label='F1', color='blue')
plt.plot(variables_exist, F2_values, 's-', label='F2', color='red')
plt.xlabel(f'{variable_name}')
plt.ylabel('F')

plt.title(f'F1 and F2 vs {variable_name}')
plt.legend()
plt.grid(True)
# ax = plt.gca()
# ax.xaxis.set_major_formatter(ticker.ScalarFormatter())
# ax.xaxis.get_major_formatter().set_scientific(False)
# ax.ticklabel_format(style='plain', axis='x')

plt.savefig(os.path.join(plot_folder, f"figure_All_F1_F2_vs_{variable_name}.png"), dpi=300)
plt.close()

# Plot uHMAX vs variable
uHMAX_values = [data[variable]["uHMAX_t"][0] for variable in data]
plt.figure(figsize=(8, 6))
plt.plot(variables_exist, uHMAX_values, '^-', label='uHMAX', color='green')
plt.xlabel(f'{variable_name}')
plt.ylabel('uHMAX')
plt.title(f'uHMAX vs {variable_name}')

plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_folder, f"figure_All_uHMAX_vs_{variable_name}.png"), dpi=300)
plt.close()

# Plot uHMAX vs variable
vRMAX_values = [data[variable]["vRMAX_t"][0] for variable in data]
plt.figure(figsize=(8, 6))
plt.plot(variables_exist, vRMAX_values, '^-', label='vRMAX', color='green')
plt.xlabel(f'{variable_name}')
plt.ylabel('vRMAX')
plt.title(f'vRMAX vs {variable_name}')

plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_folder, f"figure_All_vRMAX_vs_{variable_name}.png"), dpi=300)
plt.close()

RmaxMAX_values = [data[variable]["RmaxMAX_t"][0] for variable in data]
plt.figure(figsize=(8, 6))
plt.plot(variables_exist, RmaxMAX_values, '^-', label='Rmax', color='green')
plt.xlabel(f'{variable_name}')
plt.ylabel('Rmax')
plt.title(f'RmaxMAX vs {variable_name}')

plt.legend()
plt.grid(True)
plt.savefig(os.path.join(plot_folder, f"figure_All_RmaxMAX_vs_{variable_name}.png"), dpi=300)
plt.close()

# Plot F1_time, F2_time, uHMAX_t vs variable
F1_times = [data[variable]["F1_t"][1] for variable in data]
F2_times = [data[variable]["F2_t"][1] for variable in data]
uHMAX_times = [data[variable]["uHMAX_t"][1] for variable in data]
vRMAX_times = [data[variable]["vRMAX_t"][1] for variable in data]
RmaxMAX_times = [data[variable]["RmaxMAX_t"][1] for variable in data]
# 创建单个坐标轴
fig, ax = plt.subplots(figsize=(8, 6))

# 在同一张图上绘制
ax.plot(variables_exist, F2_times, 's-', label='F2 time', color='red')
ax.plot(variables_exist, uHMAX_times, '^-', label='uHMAX time', color='green')
ax.plot(variables_exist, F1_times, 'o-', label='F1 time', color='blue')
ax.plot(variables_exist, vRMAX_times, 'd-', label='vRMAX time', color='orange')
ax.plot(variables_exist, RmaxMAX_times, 'x-', label='RmaxMAX time', color='purple')

# 坐标轴和图例设置
ax.set_xlabel(f'{variable_name}')
ax.set_ylabel('Time')
ax.grid(True)
ax.legend(loc='best')

# 如果需要特定显示范围，可以自行打开以下语句:
# ax.set_ylim(0.3, 3.7)

# 保存并关闭
plt.savefig(os.path.join(plot_folder, f"figure_All_times_vs_{variable_name}.png"), dpi=300)
plt.close()