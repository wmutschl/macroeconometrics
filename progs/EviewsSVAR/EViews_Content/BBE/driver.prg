%s = @wfpath

cd {%s}

pageselect temp

exec {%s}bbe_f1
exec {%s}bbe_f2

close @objects

show bbe_figure_1
show bbe_figure_2
