#!/usr/bin/python
# Load the required libraries
import pandas as pd
import seaborn as sns
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# To enable writing to a file (save plot)
matplotlib.use('Agg')

def plotFig(df1, df2, df3, axes, fig, box_handles, figname ):
	df1.boxplot(column='time', by='p', showfliers=False, positions=range(df2.p.unique().shape[0]), ax=axes)
	sns.pointplot(x='p', y='time', data=df1.groupby('p', as_index=False).median(), positions=range(df2.p.unique().shape[0]), ax=axes, color='steelblue')

	df2.boxplot(column='time', by='p', showfliers=False, positions=range(df2.p.unique().shape[0]), ax=axes)
	sns.pointplot(x='p', y='time', data=df2.groupby('p', as_index=False).median(), ax=axes, color='coral')

	df3.boxplot(column='time', by='p', showfliers=False, positions=range(df3.p.unique().shape[0]), ax=axes)
	sns.pointplot(x='p', y='time', data=df3.groupby('p', as_index=False).median(), ax=axes, color='g')

	axes.legend(handles=box_handles, fontsize='small', loc=1, labelspacing=0.2)
	fig.add_axes(axes)
	axes.set_title('Scaling Factor')
	fig.add_axes(axes)
	fig.suptitle('')
	axes.set_xlabel('No. of Processes')
	fig.add_axes(axes)
	axes.set_ylabel('Time (s)')
	fig.add_axes(axes)
	fig.savefig(figname)

def plotLine(df1, df2, df3, axes, fig, line_handles, figname):
	d1=df1.groupby('p', as_index=False).median()
	d1.plot(kind='line', x='p', y='time', ax=axes)
	d2=df2.groupby('p', as_index=False).median()
	d2.plot(kind='line', x='p', y='time', ax=axes)
	d3=df3.groupby('p', as_index=False).median()
	d3.plot(kind='line', x='p', y='time', ax=axes)

	axes.legend(handles=line_handles, fontsize='small', loc=1, labelspacing=0.2)
	fig.add_axes(axes)
	axes.set(xlim=(0, 17))

	axes.set_title('Scaling Factor')
	fig.suptitle('')
	axes.set_xlabel('No. of Processes')
	axes.set_ylabel('Time (s)')
	fig.add_axes(axes)
	fig.savefig(figname)

#CSE data1
path1 = "./cse/data1"
df1 = pd.read_csv(path1+"/tt.csv")
df2 = pd.read_csv(path1+"/pt.csv")
df3 = pd.read_csv(path1+"/ppt.csv")

#CSE data2
path2 = "./cse/data2"
df4 = pd.read_csv(path2+"/tt.csv")
df5 = pd.read_csv(path2+"/pt.csv")
df6 = pd.read_csv(path2+"/ppt.csv")

#HPC data1
path3 = "./hpc/data1"
df7 = pd.read_csv(path3+"/tt.csv")
df8 = pd.read_csv(path3+"/pt.csv")
df9 = pd.read_csv(path3+"/ppt.csv")

#HPC data2
path4 = "./hpc/data2"
df10 = pd.read_csv(path4+"/tt.csv")
df11 = pd.read_csv(path4+"/pt.csv")
df12 = pd.read_csv(path4+"/ppt.csv")

# Get an axes handle from matplotlib for each of 8 figures
fig1, ax1 = plt.subplots()
fig2, ax2 = plt.subplots()
fig3, ax3 = plt.subplots()
fig4, ax4 = plt.subplots()
fig5, ax5 = plt.subplots()
fig6, ax6 = plt.subplots()
fig7, ax7 = plt.subplots()
fig8, ax8 = plt.subplots()

steelblue_box = mlines.Line2D([], [], color='steelblue', marker='o', markersize=3, label="Total Time")
coral_box = mlines.Line2D([], [], color='coral', marker='o', markersize=3, label="Avg Processing Time")
g_box = mlines.Line2D([], [], color='g', marker='o', markersize=3, label="Avg Pre-processing Time")

steelblue_line = mlines.Line2D([], [], color='steelblue', marker='_', markersize=1, label="Total Time")
coral_line = mlines.Line2D([], [], color='coral', marker='_', markersize=1, label="Avg Processing Time")
g_line = mlines.Line2D([], [], color='g', marker='_', markersize=1, label="Avg Pre-processing Time")

box_handles = [steelblue_box, coral_box, g_box ]
line_handles = [steelblue_line, coral_line, g_line ]

plotFig ( df1, df2, df3, ax1, fig1, box_handles, path1+"/plot.png" )
plotFig ( df4, df5, df6, ax2, fig2, box_handles, path2+"/plot.png" )
plotFig ( df7, df8, df9, ax3, fig3, box_handles, path3+"/plot.png" )
plotFig ( df10, df11, df12, ax4, fig4, box_handles, path4+"/plot.png" )

plotLine ( df1, df2, df3, ax5, fig5, line_handles, path1+"/plotL.png" )
plotLine ( df4, df5, df6, ax6, fig6, line_handles, path2+"/plotL.png" )
plotLine ( df7, df8, df9, ax7, fig7, line_handles, path3+"/plotL.png" )
plotLine ( df10, df11, df12, ax8, fig8, line_handles, path4+"/plotL.png" )

plt.title('Scaling Factor')
plt.suptitle('')
#plt.show()

