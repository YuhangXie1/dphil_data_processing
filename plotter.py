import matplotlib.pyplot as plt


@dataclass
class PlotConfig:
    markerstyle_map: Dict[str,str]
    markercolor_map: Dict[str,str]
    linestyle_map: Dict[str,str]
    line_color_map: Dict[Any,tuple]
    line_alpha_map: Dict[Any,float]
    plot_exclude: Dict[str, list]
    alpha_used: bool
    save_file_path: str

# Functions
def default_plot_config(save_file_path: str = ".") -> PlotConfig:
    return PlotConfig(
        markerstyle_map = {"JBL001":"s", "JBL137":"o", "media":"^"},
        markercolor_map = {"JBL001":"black", "JBL137":"blue", "media":"black"},
        linestyle_map = {"WM-met+": "solid", "WM-met-": "dashed", "LB": "solid"},
        line_color_map = {},  # build from data later
        line_alpha_map = {},
        plot_exclude = {"cells":[], "media":[], "green_intensity":[], "red_intensity":[]},
        alpha_used = False,
        save_file_path = save_file_path
    )
DefaultConfig = default_plot_config()


def plot_x_y(
        xdata : list,
        ydata : list,
        xlabel : str | None = None,
        ylabel : str | None = None,
        title : str | None = None,
        title_extra : str | None = None,
        config : PlotConfig = DefaultConfig,
        save_image : bool = False
        ):
    
    # Checking and setting defaults
    if ylabel is None:
        ylabel = ydata

    if title is None:
        title = f"{ylabel} - {plot_type} - {title_extra}"
    else:
        title = title + " - " + title_extra


    fig, axs = plt.subplots()
    for row_index, row in dataframe.iterrows():
        index_name_cells = row["cells"]
        index_name_media = row["media"]
        index_name_green_intensity = row["green_intensity"]
        index_name_red_intensity = row["red_intensity"]

        if (index_name_cells in config.plot_exclude["cells"]
            or index_name_media in config.plot_exclude["media"]
            or index_name_green_intensity in config.plot_exclude["green_intensity"]
            or index_name_red_intensity in config.plot_exclude["red_intensity"]):
            continue

        if config.alpha_used is True:
            alpha = config.line_alpha_map[index_name_green_intensity]
        else:
            alpha = 1

        if plot_type == "average":
            axs.errorbar(row["timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,

                        color = config.line_color_map[index_name_green_intensity],
                        marker = config.markerstyle_map[index_name_cells],
                        markerfacecolor = config.markercolor_map[index_name_cells],
                        markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                        linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                        alpha = alpha,
                        )
            
        elif plot_type == "all":
            for repeat in row[f"{y_data}_raw_array"]:
                axs.plot(row["timepoints"], repeat,
                        
                            color = config.line_color_map[index_name_green_intensity],
                            marker = config.markerstyle_map[index_name_cells],
                            markerfacecolor = config.markercolor_map[index_name_cells],
                            markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                            linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                            alpha = alpha,
                            )
        
        elif plot_type == "both":
            axs.errorbar(row["timepoints"], row[f"{y_data}_average"],
                        yerr = row[f"{y_data}_std"], capsize = 2.0,

                        #color = config.line_color_map[index_name_green_intensity],
                        marker = config.markerstyle_map[index_name_cells],
                        markerfacecolor = config.markercolor_map[index_name_cells],
                        markeredgecolor = config.markercolor_map[index_name_cells], markersize = 3.0,
                        linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                        alpha = 1.0,
                        )
            
            for repeat in row[f"{y_data}_raw_array"]:
                axs.plot(row["timepoints"], repeat,
                        
                            #color = config.line_color_map[index_name_green_intensity],
                            linestyle = config.linestyle_map[index_name_media], linewidth = 1.0,
                            alpha = 0.2,
                            )
                            

    leg1 = axs.legend(
        handles=cell_handles,
        title="Cell type",
        loc="upper left",
        bbox_to_anchor=(1.02, 1.00),
        borderaxespad=0.0,
    )
    axs.add_artist(leg1)

    leg2 = axs.legend(
        handles=media_handles,
        title="Media",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.65),
        borderaxespad=0.0,
    )
    axs.add_artist(leg2)

    leg3 = axs.legend(
        handles=intensity_handles,
        title="Green light intensity",
        loc="upper left",
        bbox_to_anchor=(1.02, 0.35),
        borderaxespad=0.0,
    )

    axs.add_artist(leg3)
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
