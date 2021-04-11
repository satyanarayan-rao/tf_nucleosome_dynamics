library(ggplot2)
library(grid)
library(yaml)
library(ggthemes)
library(scales)
library(stringr)
library(plotrix)
args = commandArgs(trailingOnly = T) 
options(bitmapType='cairo')

dt = read.table(args[1], sep = "", header = F, stringsAsFactor = F)
names(dt) = c("chrom_details", "cnt")
dt$cnt = as.numeric(dt$cnt)
positive = c("ENCODE", "Carroll", "positive")
negative = c("neg_ctl", "negative")
get_class = function (s) {
  s1 =  str_replace(s, "`neg_ctl", "%negative")
  return (tail(unlist(strsplit(s1, split = "%")),1))}
assign_pos_neg = function (x) {if (x %in% positive) {return ("Pos")} else if (x %in% negative) { return ("Neg")} } 
cl_and_source = apply(dt[, "chrom_details", drop = F], 1, get_class)
dt["cls"] = cl_and_source 

pos_sub_df = dt[dt$cls == "positive", ] 
# get all the positions that are > 50 
gt_50 = which(pos_sub_df$cnt > 50) 
pos_sub_df[gt_50, "cnt"] = 50
pd_pos = hist(pos_sub_df$cnt, breaks = seq(0, 50, 1), plot = F)
qq_pos = data.frame(x_tics = pd_pos$mids,  
            y_val = cumsum(pd_pos$density)/max(cumsum(pd_pos$density)), 
            cls = rep("positive", length(pd_pos$mids)) )
#plot (x = qq$x_tics, y =qq$y_val, lty = 2, col = "red", type = "l")
# fit a spline and then predict the value at 0.5
pos.spl = with(qq_pos, smooth.spline(y_val, x_tics, df = 10))
at_point_five_pos = predict(pos.spl, 0.5) 
pos.spl_inv = with(qq_pos, smooth.spline(x_tics, y_val))

xtics = seq(0, 50, 5)
xtic_labels = as.character(xtics)
xtic_labels[length(xtic_labels)] = paste0("\u2265 ", xtics[length(xtics)])
#plot (y_val ~ x_tics, data = qq_pos, main = args[4], cex = 0, xaxt = "n")
#lines(pos.spl_inv, col = "red")
#axis(side = 1, at = xtics, labels = xtic_labels)

neg_sub_df = dt[dt$cls == "negative", ] 
# get all the positions that are > 50 
gt_50 = which(neg_sub_df$cnt > 50) 
neg_sub_df[gt_50, "cnt"] = 50
pd_neg = hist(neg_sub_df$cnt, breaks = seq(0, 50, 1), plot = F)
qq_neg = data.frame(x_tics = pd_neg$mids,  
            y_val = cumsum(pd_neg$density)/max(cumsum(pd_neg$density)), 
            cls = rep("negative", length(pd_neg$mids)) ) 
neg.spl = with(qq_neg, smooth.spline(y_val, x_tics, df = 10))
at_point_five_neg = predict (neg.spl, 0.5) 

neg.spl_inv = with(qq_neg, smooth.spline(x_tics, y_val))
#par(new = TRUE)
#plot (y_val ~ x_tics, data = qq_neg, xlab = "", ylab = "", axes = FALSE, cex = 0)
#lines(neg.spl_inv, col = "blue")

d = round(at_point_five_pos$y - at_point_five_neg$y, 2)
#text (x = 30, y = 0.1, paste0("Diff at 0.5 = ", d), pos = 4)
#segments(at_point_five_neg$y, 0.5, at_point_five_pos$y, 0.5, lty =2, lwd = 2)
#abline(h = 0.5, col = "#636363", lty = 2)

#legend(30, 0.4, legend = c("Positive Sites", "Negative Sites"), 
#      col = c("Red", "Blue"), lty = c(1,1), cex = 0.8)


to_plot_df = rbind (qq_pos, qq_neg)

to_save_df = to_plot_df
to_save_df$system = args[4]
write.table(to_save_df, args[6], row.names = F, 
            col.names = T, quote = F, sep = "\t")
png(args[3], width = 6, height = 5, units = "in", res = 150)

plt = ggplot (to_plot_df, aes(x = x_tics, y = y_val, colour = cls)) + geom_line() + geom_rangeframe() + theme_few()
plt + scale_x_continuous(breaks = seq(0,50,5), labels = xtic_labels) + 
  xlab("Number of fragments in TFBS") + ylab("eCDF") + 
  ggtitle(args[4]) + 
  theme(plot.title = element_text(hjust = 0.5)) + 
  geom_hline(yintercept = 0.5, lty = 2, col = "#636363") + 
  annotate("text", x = 40, y = 0.4, label = paste0("Difference at 0.5 = ", d)) +
  scale_colour_manual(breaks = c("positive", "negative"), 
                      values = c("red", "blue"), labels = c("Positive", "Negative")) + 
  theme(legend.title = element_blank(), legend.position = c(0.8, 0.2)) + 
  ylim(c(0,1)) 
dev.off ()

write.table(dt, args[5], row.names = F, 
            col.names = T, quote = F, sep = "\t")

png (args[2], height = 5, width = 6, units = "in", res = 150)

plt = ggplot(dt, aes (x = cnt, fill = cls)) + 
      geom_histogram(color="#e9ecef", alpha=0.5, position = 'identity', bins = 20) + 
      scale_fill_manual(values=c("#69b3a2", "#404080")) + geom_rangeframe() + 
      theme_few() + ggtitle (args[4]) + theme(plot.title = element_text(hjust=0.5)) + 
      xlab("Number of fragments in TFBS") + ylab("Count") + xlim(c(0,50)) 
print (plt)
dev.off()

