from Polarization import *

r"""
Grouping of rcParams for matplotlib. For easily changing global settings such as ticks, legend placement, lenght and fontsize.
As it currently stands there is no reason to put in a separate params file.
"""
matplotlib.rcParams['figure.figsize'] = 12, 8
matplotlib.rcParams['font.size']= 18
matplotlib.rcParams['legend.handlelength']= 4
matplotlib.rcParams['legend.loc'] = "best"
matplotlib.rcParams['figure.subplot.hspace'] = 0.05
matplotlib.rcParams["xtick.minor.visible"] =  True
matplotlib.rcParams["ytick.minor.visible"] =  True
matplotlib.rcParams["lines.linewidth"] =  4
matplotlib.rcParams["xtick.top"] =  True
matplotlib.rcParams["ytick.right"] =  True
#matplotlib.pyplot.annotate(r"$q=\frac{\theta}{\theta_j}$", (0.91, 0.9), xycoords='figure fraction') #(Annotations on the plot e.g. equations)
#.legend(loc='center left', bbox_to_anchor=(1, 0.5)) Put legend box outside of the graph


Folder_dicts = {"folders":[target_folder + "Output/Ordered/", target_folder + "Output/Perpendicular/",
						   target_folder + "Output/Parallel/", target_folder + "Output/Toroidal/"],
 				"plot_folders":[target_folder + "Output/Ordered/Plots/", target_folder + "Output/Perpendicular/Plots/",
								target_folder + "Output/Parallel/Plots/", target_folder + "Output/Toroidal/Plots/"]}


r"""
Read the output files obtained from the Polarization.py code and return as a pandas dictionary for plotting
"""
Ordered_Output = pandas.read_csv(Folder_dicts["folders"][0] + "Ordered_Output.csv", sep = "\t", index_col = 0)
Ordered_Output_2 = pandas.read_csv("Batch_2/Output/Ordered/" + "Ordered_Output.csv", sep = "\t", index_col = 0)
Perpendicular_Output = pandas.read_csv(Folder_dicts["folders"][1] + "Perpendicular_Output.csv", sep= "\t", index_col = 0)
Parallel_Output = pandas.read_csv(Folder_dicts["folders"][2] + "Parallel_Output.csv", sep = "\t", index_col = 0)
Toroidal_Output = pandas.read_csv(Folder_dicts["folders"][3] + "Toroidal_Output.csv", sep = "\t", index_col = 0)
t_dec = pandas.read_csv(target_folder + "t_dec_output.csv", sep = "\t", index_col = 0)

Ordered_Output_sum = Ordered_Output.groupby([r'$t$']).sum()
print(Ordered_Output)
Ordered_Output_2_sum = Ordered_Output_2.groupby([r'$t$']).sum()
print(Ordered_Output_2)
absd = Ordered_Output_sum + Ordered_Output_2_sum
absd = (absd-absd.min())/(absd.max()-absd.min())
absd = absd.reset_index()
print(absd)

r"""
Plotting of the dictionaries. 
Not included in rcParams:
.legend(loc = ,bbox_to_anchor =); extra function from pandas that allows to put the legend box outside of the graph. Useful for multiple runs of q.
.xscale('log'); logarithmic scale for the x axis. Important as the time distribution is currently a geometric progression (log) so we can explore the early points.
"""
plot_ord =	Ordered_Output.plot(x=r'$t$',
	style = ['--', '-.', '-', ':'],
#	yticks = numpy.linspace(-0.2, 0.5, 7),
#	xticks = numpy.linspace(0.0, 3, 7)
	).legend(loc='center left', bbox_to_anchor=(1, 0.5))

for k in range(t_dec[r'$t_{{dec}}$'].shape[0]):
	#this is for adding vertical lines (corresponding to the deceleration timescale) to the plot. Will change to a single line when working only with 1 value of angle instead of a list of angles
	matplotlib.pyplot.axvline(t_dec[r'$t_{{dec}}$'][k], linewidth=2, color="black")

matplotlib.pyplot.xscale('log')
fig_ord = plot_ord.get_figure()
fig_ord.savefig(Folder_dicts["plot_folders"][0] + "Output_Ordered.png", bbox_inches = "tight")

plot_ord =	absd.plot(x=r'$t$',
	style = ['--', '-.', '-', ':'],
#	yticks = numpy.linspace(-0.2, 0.5, 7),
#	xticks = numpy.linspace(0.0, 3, 7)
	).legend(loc='center left', bbox_to_anchor=(1, 0.5))

for k in range(t_dec[r'$t_{{dec}}$'].shape[0]):
	#this is for adding vertical lines (corresponding to the deceleration timescale) to the plot. Will change to a single line when working only with 1 value of angle instead of a list of angles
	matplotlib.pyplot.axvline(t_dec[r'$t_{{dec}}$'][k], linewidth=2, color="black")

matplotlib.pyplot.xscale('log')
fig_ord = plot_ord.get_figure()
fig_ord.savefig(Folder_dicts["plot_folders"][0] + "Output_Ordered_sum.png", bbox_inches = "tight")




plot_per = Perpendicular_Output.plot(x=r'$t$',
	style = ['--', '-.', '-', ':'],
#	yticks = numpy.linspace(-0.2, 0.5, 7),
#	xticks = numpy.linspace(0.0, 3, 7)
	).legend(loc='center left', bbox_to_anchor=(1, 0.5))

for k in range(t_dec[r'$t_{{dec}}$'].shape[0]):
	#this is for adding vertical lines (corresponding to the deceleration timescale) to the plot. Will change to a single line when working only with 1 value of angle instead of a list of angles
	matplotlib.pyplot.axvline(t_dec[r'$t_{{dec}}$'][k], linewidth=2, color="black")

matplotlib.pyplot.xscale('log')
fig_per = plot_per.get_figure()
fig_per.savefig(Folder_dicts["plot_folders"][1] + "Output_Perpendicular.png", bbox_inches = "tight")


plot_par = Parallel_Output.plot(x=r'$t$',
	style = ['--', '-.', '-', ':'],
#	yticks = numpy.linspace(-0.2, 1.0, 7),
#	xticks = numpy.linspace(0.0, 3, 7)
	).legend(loc='center left', bbox_to_anchor=(1, 0.5))

for k in range(t_dec[r'$t_{{dec}}$'].shape[0]):
	#this is for adding vertical lines (corresponding to the deceleration timescale) to the plot. Will change to a single line when working only with 1 value of angle instead of a list of angles
	matplotlib.pyplot.axvline(t_dec[r'$t_{{dec}}$'][k], linewidth=2, color="black")

matplotlib.pyplot.xscale('log')
fig_par = plot_par.get_figure()
fig_par.savefig(Folder_dicts["plot_folders"][2] + "Output_Parallel.png", bbox_inches = "tight")



plot_tor = Toroidal_Output.plot(x=r'$t$',
	style = ['--', '-.', '-', ':'],
#	yticks = numpy.linspace(0.0, 1.0, 7),
	#xticks = numpy.linspace(0.0, 3, 7)
	).legend(loc='center left', bbox_to_anchor=(1, 0.5))

for k in range(t_dec[r'$t_{{dec}}$'].shape[0]):
	#this is for adding vertical lines (corresponding to the deceleration timescale) to the plot. Will change to a single line when working only with 1 value of angle instead of a list of angles
	matplotlib.pyplot.axvline(t_dec[r'$t_{{dec}}$'][k], linewidth=2, color="black")

matplotlib.pyplot.xscale('log')
fig_tor = plot_tor.get_figure()
fig_tor.savefig(Folder_dicts["plot_folders"][3] + "Output_Toroidal.png", bbox_inches = "tight")


