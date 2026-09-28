from matplotlib.lines import Line2D
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from typing import Dict, Any
from dataclasses import dataclass
import os
from os import path
from pathlib import Path

#===================================================================================#
#Headers are LED,FOV,Frame,Mean_Brightness,Median_Brightness

filepath = '26-09-25 frame brightness mm/frame_brightness.csv'
light_regime_auto = True
light_regime_custom = {"green" : [0,1,2,3,4,5],
                    "dark" : [6,7,8,9,10,11],
                    "red" : [12,13,14,15,16,17],}



def main(filepath):
    dataframe = pd.read_csv(filepath, header = 0, index_col= None)
    
    plot_exclude = {"fov": [19],
                    "light_regime": []}

    
    DefaultConfig.plot_exclude = plot_exclude
    plot_brightness_by_frame(dataframe, "LED450NM")
    


#===================================================================================#
@dataclass
class PlotConfig:
    markerstyle_map: Dict[str,str]
    markercolor_map: Dict[str,str]
    linestyle_map: Dict[str,str]
    line_color_map: Dict[Any,tuple]
    line_alpha_map: Dict[Any,float]
    plot_exclude: Dict[str, list]
    alpha_used: bool
    light_regime_auto: bool
    light_regime: Dict[str, list]
    save_file_path: str

    default_markerstyle: str = "o"
    default_markercolor: str = "blue"
    default_linestyle: str = "solid"
    default_linecolor: str = "black"
    default_linealpha: float = 1.0

def default_plot_config(save_file_path: str = ".") -> PlotConfig:
    return PlotConfig(
        markerstyle_map = {},
        markercolor_map = {},
        linestyle_map = {},
        line_color_map = {"green": (0.0, 1.0, 0.0), "dark": (0.5, 0.5, 0.5), "red": (1.0, 0.0, 0.0)},
        line_alpha_map = {},
        plot_exclude = {"fov": [], "light_regime": []},
        alpha_used = False,
        light_regime_auto = True,
        light_regime = {},
        save_file_path = save_file_path
    )
DefaultConfig = default_plot_config()

def add_light_regime_to_dataframe(dataframe, light_regime_auto = True, light_regime_custom = None):
    fov_count = len(dataframe["FOV"].unique())
    fov_sections = 3
    fov_sections_split = [1,1,1] #relative size of sections [Green, Dark, Red]
    number_of_FOVs_per_section = int(fov_count / fov_sections) #middle section is dark, otherwise green and red
    green_regime = []
    dark_regime = []
    red_regime = []
    other_regime = []
    
    for fov in dataframe["FOV"].unique():
        if number_of_FOVs_per_section is not None:
            if fov < (number_of_FOVs_per_section * fov_sections_split[0]):
                green_regime.append(fov)
            elif fov >= (number_of_FOVs_per_section * fov_sections_split[0]) and fov < (number_of_FOVs_per_section * (fov_sections_split[0]+fov_sections_split[1])):
                dark_regime.append(fov)
            elif fov >= (number_of_FOVs_per_section * (fov_sections_split[0]+fov_sections_split[1])) and fov < (number_of_FOVs_per_section * (fov_sections_split[0]+fov_sections_split[1]+fov_sections_split[2])):
                red_regime.append(fov)
            else:
                other_regime.append(fov)

    print(f"Max FOVs: {fov_count}, Green: {list(green_regime)}, Dark: {list(dark_regime)}, Red: {list(red_regime)}, None: {list(other_regime)}")

    if light_regime_auto:
        dataframe["light_regime"] = "dark"
        dataframe.loc[dataframe["FOV"].isin(green_regime), "light_regime"] = "green"
        dataframe.loc[dataframe["FOV"].isin(dark_regime), "light_regime"] = "dark"
        dataframe.loc[dataframe["FOV"].isin(red_regime), "light_regime"] = "red"
    else:
        dataframe["light_regime"] = "dark"
        for regime, fovs in light_regime_custom.items():
            dataframe.loc[dataframe["FOV"].isin(fovs), "light_regime"] = regime
    return dataframe


def plot_diff_by_fov(dataframe, LED, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Frames",
                            config: PlotConfig = DefaultConfig, save_image = False,
                           ):
    """
    Plots brightness over frames
    
    :param dataframe: dataframe where the data is (raw)
    :param LED: LED to plot (e.g. 450nm)
    :param ylabel: optional override for the ylabel text
        :type ylabel: str | None
    :param title: optional override for title text
        :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
        :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
        :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """

    if ylabel is None:
        ylabel = "Brightness (AU)"

    if title is None:
        title = f"Average pixel brightness by FOV - {title_extra}"
    else:
        title = title + " - " + title_extra

    #slicing data
    dataframe = dataframe.loc[dataframe["LED"] == LED].copy()
    dataframe = add_light_regime_to_dataframe(dataframe)

    fig, axs = plt.subplots()
    for fov, group in dataframe.groupby("FOV"):
        light_regime = group["light_regime"].iloc[0]
        if (fov in config.plot_exclude["fov"] or light_regime in config.plot_exclude["light_regime"]):
            continue

        #calculating difference in brightness
        bright_min = group["Mean_Brightness"].min()
        bright_max = group["Mean_Brightness"].max()
        bright_end = group["Mean_Brightness"][-1]
        bright_end_period = 10 #frames to average over
        bright_end_average = np.mean(group["Mean_Brightness"][-bright_end_period:-1])
        
        x = fov
        y = (bright_max - bright_min, bright_end - bright_min, bright_end_average - bright_min)

        axs.grouped_bar(y, tick_label = x, color = config.line_color_map[light_regime],
                    #linestyle = config.linestyle_map[fov],
                    linewidth = 1.0,
                    )
        
        light_regime_handle = [
                Line2D(
                    [], [], 
                    color = config.line_color_map[light_regime],
                    label = light_regime
                )
                for light_regime in config.line_color_map
            ]

        leg1 = axs.legend(
                handles=light_regime_handle,
                title="Light regime",
                loc="upper left",
                bbox_to_anchor=(1.02, 1.00),
                borderaxespad=0.0,
            )
        axs.add_artist(leg1)

        # FOV label at start
        axs.annotate(
            str(fov),
            xy=(x.iloc[0], y.iloc[0]),
            xytext=(-5, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize = 6,
        )

        # FOV label at end
        axs.annotate(
            str(fov),
            xy=(x.iloc[-1], y.iloc[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize = 6,
        )

    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout(rect=[0, 0, 0.75, 1])

    if save_image is True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)

def plot_brightness_by_frame(dataframe, LED, ylabel: str | None = None, title: str | None = None, title_extra: str = "",
                            xlabel = "Frames",
                            config: PlotConfig = DefaultConfig, save_image = False,
                           ):
    """
    Plots brightness over frames
    
    :param dataframe: dataframe where the data is (raw)
    :param LED: LED to plot (e.g. 450nm)
    :param ylabel: optional override for the ylabel text
        :type ylabel: str | None
    :param title: optional override for title text
        :type title: str | None
    :param title_extra: optional extra to tag onto the end of the title text
        :type title_extra: str
    :param xlabel: optional override for the xlabel text
    :param config: config file for styling and saving the image
        :type config: PlotConfig
    :param save_image: if True, then image is saved instead of displayed
    """

    if ylabel is None:
        ylabel = "Brightness (AU)"

    if title is None:
        title = f"Average pixel brightness over frames - {title_extra}"
    else:
        title = title + " - " + title_extra

    #slicing data
    dataframe = dataframe.loc[dataframe["LED"] == LED].copy()
    dataframe = add_light_regime_to_dataframe(dataframe)

    fig, axs = plt.subplots()
    for fov, group in dataframe.groupby("FOV"):
        light_regime = group["light_regime"].iloc[0]
        if (fov in config.plot_exclude["fov"] or light_regime in config.plot_exclude["light_regime"]):
            continue
        
        x = group["Frame"]
        y = group["Mean_Brightness"]

        axs.plot(x, y, color = config.line_color_map[light_regime],
                    #linestyle = config.linestyle_map[fov],
                    linewidth = 1.0,
                    )
        
        light_regime_handle = [
                Line2D(
                    [], [], 
                    color = config.line_color_map[light_regime],
                    label = light_regime
                )
                for light_regime in config.line_color_map
            ]

        leg1 = axs.legend(
                handles=light_regime_handle,
                title="Light regime",
                loc="upper left",
                bbox_to_anchor=(1.02, 1.00),
                borderaxespad=0.0,
            )
        axs.add_artist(leg1)

        # FOV label at start
        axs.annotate(
            str(fov),
            xy=(x.iloc[0], y.iloc[0]),
            xytext=(-5, 0),
            textcoords="offset points",
            ha="right",
            va="center",
            fontsize = 6,
        )

        # FOV label at end
        axs.annotate(
            str(fov),
            xy=(x.iloc[-1], y.iloc[-1]),
            xytext=(5, 0),
            textcoords="offset points",
            ha="left",
            va="center",
            fontsize = 6,
        )

    axs.set_xlabel(xlabel)
    axs.set_ylabel(ylabel)
    axs.set_title(title)
    fig.tight_layout(rect=[0, 0, 0.75, 1])

    if save_image is True:
        try:
            Path(os.path.join(config.save_file_path, "figures")).mkdir(parents = True, exist_ok = True)
            save_title = title.replace("/","_div_")
            plt.savefig(os.path.join(config.save_file_path, "figures", f"{save_title}.png"))
            plt.close(fig)
        except:
            print("Could not save image, filepath not valid")
    else:
        plt.show()
        plt.close(fig)

if __name__ == "__main__":
    main(filepath)