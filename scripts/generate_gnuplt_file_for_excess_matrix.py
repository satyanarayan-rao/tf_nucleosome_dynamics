import os
import sys
inp_fp = open(sys.argv[1])
base_plt = open (sys.argv[2])
out_fp = open(sys.argv[3], "w")
total_rows = int(sys.argv[4]) - 1
out_fname = sys.argv[5]

significant_str_fp = open(sys.argv[8])
for line in base_plt:
    out_fp.write(line)

base_plt.close()
cnt = 0
cell_cnt = 0
max_val = -1000 
min_val = 1000
for line, s_l in zip(inp_fp,significant_str_fp):
    l_items = line.strip().split()
    v_cnt = 0
    s_items = s_l.strip().split()
    for v, s in zip(l_items, s_items):
        if s !="n.s.":
            #to_write = "set label %s at %s, %s '%s^{ %s }' center front"%(cell_cnt+1, v_cnt, total_rows - cnt, v, s)
            to_write = "set label %s at %s, %s '%s' center front"%(cell_cnt+1, v_cnt, total_rows - cnt, s)
            
        #out_fp.write("set label {lab_cnt} at {x}, {y} '{v}^{{{}}{s}{{}}}' center front\n".\
        #             format(lab_cnt = cell_cnt+1, x = v_cnt, y = total_rows - cnt,
        #                    v = v, s = s) ) 
        out_fp.write(to_write + "\n")
        v_cnt +=1
        cell_cnt +=1
        f_v = float(v)
        if f_v > max_val: 
            max_val = f_v 
        if f_v < min_val:
            min_val = f_v 
    cnt +=1
#out_fp.write("set cbra [{min_v}:{max_v}]\n".format(min_v=min_val, max_v=max_val)) 
out_fp.write("set palette defined (0.5 '#1f78b4', 1 '#ffffff', 1.5 '#ef8a62')\n") 
out_fp.write("set cbra [0.5:1.5]\n".format(min_v=min_val, max_v=max_val)) 
out_fp.write("set cbtics ('0.5' 0.5, '1.0' 1, '1.5' 1.5)\n")

out_fp.write("set output '{o}'\n".format(o=out_fname))

out_fp.write("set yra [-0.5:{v}]\n".format(v=total_rows+0.5))
out_fp.write("set xra [-0.5:{v}]\n".format(v=total_rows+0.5))

ylab_list = ["\"" + "Cl" + str(i + 1) + "\"" + " " + str(total_rows - i) for i in range(total_rows, -1, -1) ]

ylab_list_str = ",".join(ylab_list)

xlab_list = ["\"" + "Cl" + str(i + 1) + "\"" + " " + str(i) for i in range(total_rows+1) ]
xlab_list_str = ",".join(xlab_list)
out_fp.write("set x2tics (" + xlab_list_str + ")\n")
out_fp.write("set noxtics\n")
out_fp.write("set ytics (" + ylab_list_str + ")\n")
out_fp.write("set ylab " + "'" + sys.argv[6] +"'" + ' noenhanced' + "\n")
out_fp.write("set x2lab " + "'" + sys.argv[7] + "'" + ' noenhanced' +  "\n")

out_fp.write("plot '<tac {f}' matrix with image notitle".format(f=sys.argv[1]))

out_fp.close()
inp_fp.close()
