import pandas as pd
import matplotlib.pyplot as plt

headers = [
    "Experiment time (s)",
    "Measured OD",
    "Turbidostat OD setpoint",
    "OD zero setpoint",
    "Thermostat setpoint (C)",
    "Heating rate",
    "Internal air temperature (C)",
    "External air temperature (C)",
    "Media temperature (C)",
    "Optogenetic actuation intensity",
    "Pump 1 rate",
    "Pump 2 rate",
    "Pump 3 rate",
    "Pump 4 rate",
    "Volume",
    "Stirring rate",
    "395nm LED set-point",
    "?1",
    "?2",
    "?3",
    "457nm LED set-point",
    "500nm LED set-point",
    "523nm LED set-point",
    "595nm LED set-point",
    "623nm LED set-point",
    "6500K LED set-point",
    "Laser set-point",
    "UV LED intensity",
    "Fluorescent protein 1 baseband",
    "Fluorescent protein 1 emission band 1",
    "Fluorescent protein 1 emission band 2",
    "Fluorescent protein 2 baseband",
    "Fluorescent protein 2 emission band 1",
    "Fluorescent protein 2 emission band 2",
    "Fluorescent protein 3 baseband",
    "Fluorescent protein 3 emission band 1",
    "Fluorescent protein 3 emission band 2",
    "Custom program parameter 1",
    "Custom program parameter 2",
    "Custom program parameter 3",
    "Custom program status",
    "Zigzag target",
    "Growth rate",
]

def plot_timecourse(chosen_data_df, columns, title = "", ylabel = "", xlabel = "", time_step = "hr", slice = None, LED_periods = None):
    
    if time_step == "hr":
        x = chosen_data_df["Experiment time (hr)"]
    else:
        x = chosen_data_df["Experiment time (s)"]

    fig, axs = plt.subplots(figsize=(10, 5))

    if LED_periods is not None:
        for start, end, color in LED_periods:
            axs.axvspan(start, end, color = color, alpha = 0.15, zorder = 0)

    for c in columns:
        axs.plot(x, chosen_data_df[c], label = c)
    
    axs.set_title(title)
    axs.set_ylabel(ylabel)
    axs.set_xlabel(xlabel)
    axs.legend(bbox_to_anchor=(1.02, 1.00))
    fig.tight_layout(rect=[0, 0, 0.85, 1])
    
    plt.show()

filepath = "26-02-05_chi_bio_JBL137_1/2026-01-30 15_46_01_M0_data.csv"
data_df = pd.read_csv(filepath, names = headers, index_col= False)
print(data_df)

data_df.insert(data_df.columns.get_loc("Experiment time (s)")+1, "Experiment time (hr)", data_df["Experiment time (s)"]/3600)

print(data_df)

#print(data_df["Experiment time (s)"])
#print(data_df["Measured OD"])

""" data_columns = headers[1:]
for column in data_columns:

    fig, axs = plt.subplots()

    axs.plot(data_df["Experiment time (s)"], data_df[column])
    axs.set_title(column)
    axs.set_ylabel(column)
    axs.set_xlabel("time (s)")

    plt.show() """



# Diagnostic figures
plot_timecourse(data_df, ["Internal air temperature (C)", "External air temperature (C)", "Media temperature (C)", "Thermostat setpoint (C)"], title= "Temperature", ylabel="Temperature (C)", xlabel = "Time (s)")


# Figures
LED_set = [(0,12,"red"),(12,24,"green"),(24,36,"red"),(36,48,"green"),(48,60,"red")]
#plot_timecourse(data_df, ["Measured OD", "Turbidostat OD setpoint"], title= "Measured OD", ylabel="OD600", xlabel = "Time (hr)", LED_periods = LED_set)
#plot_timecourse(data_df, ["Fluorescent protein 1 emission band 1"], title= "Fluorescent protein 1 emission band 1", ylabel="GFP395", xlabel = "Time (hr)", LED_periods = LED_set)
#plot_timecourse(data_df, ["Fluorescent protein 1 emission band 2"], title= "Fluorescent protein 1 emission band 2", ylabel="GFP395", xlabel = "Time (hr)", LED_periods = LED_set)
#plot_timecourse(data_df, ["Fluorescent protein 1 baseband"], title= "Fluorescent protein 1 baseband", ylabel="GFP395", xlabel = "Time (hr)", LED_periods = LED_set)
#plot_timecourse(data_df, ["Growth rate"], title= "Growth rate", ylabel="Growth rate", xlabel = "Time (hr)", LED_periods = LED_set)


# Custom figures

# Growth rate overlaid with growth curves and fluorescence
fig, axs = plt.subplots(figsize=(10, 5))
axs2 = axs.twinx()
axs3 = axs.twinx()
axs.spines["left"].set_color("green")
axs2.spines["right"].set_color("blue")
axs3.spines["right"].set_position(("outward", 60))
axs3.spines["right"].set_color("red")

for start, end, color in LED_set:
    axs.axvspan(start, end, color = color, alpha = 0.15, zorder = 0)

axs.plot(data_df["Experiment time (hr)"], data_df["Fluorescent protein 1 emission band 1"], label = "Fluorescence", color = "green")
axs.set_ylabel("Fluorescence", color = "green")

axs2.plot(data_df["Experiment time (hr)"], data_df["Measured OD"], label = "Measured OD", color = "blue")
axs2.set_ylabel("Measured OD", color = "blue")

axs3.plot(data_df["Experiment time (hr)"], data_df["Growth rate"], label = "Growth rate", color = "red")
axs3.set_ylabel("Growth rate", color = "red")

axs.set_title("Growth rate, OD, fluorescence")
axs.set_xlabel("Time (hr)")

lines1, labels1 = axs.get_legend_handles_labels()
lines2, labels2 = axs2.get_legend_handles_labels()
lines3, labels3 = axs3.get_legend_handles_labels()
axs.legend(
    lines1 + lines2 + lines3,
    labels1 + labels2 + labels3,
    bbox_to_anchor=(1.02, 1.20)
)

fig.tight_layout(rect=[0, 0, 0.85, 1])

plt.show()

