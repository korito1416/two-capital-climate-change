
bash ./conduction/baseline_plot.sh

## Figure 6,7 and Table 2 FK
bash ./conduction/FeymannKacs_plot_composite.sh


# Table 3 should change to grid 0.1
bash ./conduction/DateValueReport_table.sh


# Figure 8
bash ./conduction/StocasticTrajectories_plot.sh "true" "false" "false"
sleep 1800
# Figure 11
bash ./conduction/StocasticTrajectories_plot.sh "false" "true" "false"
sleep 1800
# Figure 12
bash ./conduction/StocasticTrajectories_plot.sh "false" "false" "true"

# Table 4
bash ./conduction/FeymannKacs_plot_composite_NewPlug.sh

# Figure 10
bash ./conduction/FeymannKacs_plot_NewPlug.sh