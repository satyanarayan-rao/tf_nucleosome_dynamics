import sys
import subprocess

tsv_file=sys.argv[1]
basic_commands_fp = open (sys.argv[2])
out_fp = open (sys.argv[3], "w")
tsv_fp = open (tsv_file)

#line_breaks_file=sys.argv[2]
out = subprocess.Popen(['wc', '-l', tsv_file], stdout=subprocess.PIPE, stderr=subprocess.STDOUT)  
stdout, stderr = out.communicate()
num_lines = int (stdout.decode().split()[0]) - 1  
line_headers = tsv_fp.readline().strip().split()
len_line_headers = len(line_headers)

# what we want is the string of this format: ("x" 0, "y" 1)
line_headers_with_quotes = [ '"' + i + '"' for i in line_headers]
gnuplot_xlab_string = [ j + " " + str(i) for i,j in zip (list(range(len_line_headers)), line_headers_with_quotes )]

xtics_str = "(" + ",".join(gnuplot_xlab_string) + ")" + " " + "rotate by -45"

for line in basic_commands_fp: 
    out_fp.write(line)
out_fp.write ("set output '{f}'\n".format(f=sys.argv[4]))
out_fp.write ("set yrange [0.5:{y}]\n".format(y = num_lines + 0.5)) 
out_fp.write("set ytics 1,1,{nlines}\n".format(nlines = num_lines))
out_fp.write ("set xtics {labels}\n".format (labels = xtics_str))
out_fp.write ("set title '{f}'\n".format(f=sys.argv[5]))
out_fp.write ("plot '{f}' matrix with image notitle\n".format(f=sys.argv[1]))

basic_commands_fp.close() 
out_fp.close()
