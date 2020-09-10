Cloning Success Rates
================
Kenneth Matreyek

``` r
input_data <- read.csv(file = "MatreyekLab_Hifi_Reactions - colonies_success_rates.csv", header = T, stringsAsFactors = F)
input_data$comb_rate <- as.numeric(input_data$comb_rate)

gibson <- input_data %>% filter(!is.na(comb_rate)) %>% group_by(gibson, parts) %>% add_tally() %>% summarize(mean = mean(comb_rate), sd = sd(comb_rate), n = median(n))
```

    ## `summarise()` regrouping output by 'gibson' (override with `.groups` argument)

``` r
parts1_gibson <- input_data %>% filter(!is.na(parts) & !is.na(tested) & parts == 1 & gibson == "yes")
parts1_gibson <- data.frame("comb_rate" = parts1_gibson$comb_rate, "parts" = 1, "gibson" = "yes")
parts2_gibson <- input_data %>% filter(!is.na(parts) & !is.na(tested) & parts == 2 & gibson == "yes")
parts2_gibson <- data.frame("comb_rate" = parts2_gibson$comb_rate, "parts" = 2, "gibson" = "yes")
parts1_nogibson1 <- input_data %>% filter(!is.na(parts) & !is.na(tested) & parts == 1 & gibson == "no") 
parts1_nogibson1[parts1_nogibson1$muts_of == "","muts_of"] <- seq(1,sum(parts1_nogibson1$muts_of == ""))
parts1_nogibson <- parts1_nogibson1 %>% group_by(muts_of) %>% summarize(comb_rate = mean(comb_rate), parts = unique(parts), gibson = unique(gibson)) %>% select(-muts_of)
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
parts1_nogibson <- data.frame("comb_rate" = parts1_nogibson$comb_rate, "parts" = 1, "gibson" = "no")
parts2_nogibson <- input_data %>% filter(!is.na(parts) & !is.na(tested) & parts == 2 & gibson == "no")
parts2_nogibson <- data.frame("comb_rate" = parts2_nogibson$comb_rate, "parts" = 2, "gibson" = "no")

combined_frame <- rbind(parts1_gibson, parts2_gibson, parts1_nogibson, parts2_nogibson)
combined_frame[combined_frame$parts == 1,"parts"] <- "1-part"
combined_frame[combined_frame$parts == 2,"parts"] <- "2-parts"

Gibson_plot <- ggplot() + theme_bw() + theme(panel.grid.major.x = element_blank(), panel.grid.minor.x = element_blank()) + 
  geom_density(data = combined_frame, aes(x = comb_rate, y = ..count.., fill = gibson), alpha = 0.3, 
position="stack") +
  facet_grid(rows = vars(parts)) + xlab("Rate of success") + ylab("Number of reactions") +
  geom_vline(xintercept = 0.75, linetype = 2, color = "blue") +
  geom_vline(xintercept = 0.25, linetype = 2, color = "red")
ggsave(file = "Gibson_plot.pdf", Gibson_plot, height = 3, width = 4)
```

    ## Warning: Removed 4 rows containing non-finite values (stat_density).

``` r
Gibson_plot
```

    ## Warning: Removed 4 rows containing non-finite values (stat_density).

![](Cloning_success_rates_files/figure-gfm/Part%201-1.png)<!-- -->

``` r
# Make another set of plots using the above data to show how often it "works immediately", "works with some struggle", and simply "does not work"

temp <- ggplot_build(ggplot() + geom_density(data = combined_frame, aes(x = comb_rate, y = ..count.., fill = gibson), alpha = 0.3))
```

    ## Warning: Removed 4 rows containing non-finite values (stat_density).

``` r
part2 <- data.frame(temp$data) %>% filter(group == 2)
part1 <- data.frame(temp$data) %>% filter(group == 1)

With_gibson_density_plot <- ggplot() + 
  theme_bw() + theme(panel.grid.minor = element_blank()) +
  scale_x_continuous(limits = c(0,1), expand = c(0,0)) +
  scale_y_continuous(limits = c(0,1.1), expand = c(0,0)) +
  labs(x = "Frequency of successful clones per Gibson attempt", y = "Relative density\n(fraction of maximum)") +
  geom_point(data = part1, aes(x = x, y = density), color = "#F8766D") #with Gibson
With_gibson_density_plot
```

![](Cloning_success_rates_files/figure-gfm/Part%20two-1.png)<!-- -->

``` r
ggsave(file = "With_gibson_density_plot.pdf", With_gibson_density_plot, height = 3, width = 5)
ggsave(file = "/Users/kmatreyek/Dropbox/Website/Plots/With_gibson_density_plot.png", With_gibson_density_plot, height = 3, width = 5)

#ggplot() + geom_point(data = part2, aes(x = x, y = density), color = "#00BFC4") #without Gibson
probs <- data.frame("probability" = sample_n(part1, 100, weight = part1$density)$x)

output_frame <- data.frame("clones" = rep(seq(0,5),50), "successes_in_hundred" = 0)
for(x in 1:nrow(output_frame)){
  result_vector <- c()
  for(y in seq(1,100)){
    temp_probability <- sample_n(part1, 1, weight = part1$density)$x
    temp_clone_results <- rbinom(output_frame$clones[x], 1, prob = temp_probability)
    temp_success <- 1 %in% temp_clone_results
    result_vector <- c(result_vector, temp_success)
  }
  output_frame[x,"successes_in_hundred"] <- sum(result_vector)
}
output_frame
```

    ##     clones successes_in_hundred
    ## 1        0                    0
    ## 2        1                   46
    ## 3        2                   64
    ## 4        3                   66
    ## 5        4                   75
    ## 6        5                   78
    ## 7        0                    0
    ## 8        1                   34
    ## 9        2                   51
    ## 10       3                   71
    ## 11       4                   70
    ## 12       5                   79
    ## 13       0                    0
    ## 14       1                   41
    ## 15       2                   59
    ## 16       3                   59
    ## 17       4                   69
    ## 18       5                   74
    ## 19       0                    0
    ## 20       1                   46
    ## 21       2                   64
    ## 22       3                   58
    ## 23       4                   73
    ## 24       5                   72
    ## 25       0                    0
    ## 26       1                   44
    ## 27       2                   55
    ## 28       3                   74
    ## 29       4                   66
    ## 30       5                   67
    ## 31       0                    0
    ## 32       1                   45
    ## 33       2                   58
    ## 34       3                   73
    ## 35       4                   70
    ## 36       5                   74
    ## 37       0                    0
    ## 38       1                   45
    ## 39       2                   54
    ## 40       3                   73
    ## 41       4                   82
    ## 42       5                   77
    ## 43       0                    0
    ## 44       1                   36
    ## 45       2                   60
    ## 46       3                   67
    ## 47       4                   75
    ## 48       5                   76
    ## 49       0                    0
    ## 50       1                   45
    ## 51       2                   54
    ## 52       3                   64
    ## 53       4                   78
    ## 54       5                   74
    ## 55       0                    0
    ## 56       1                   40
    ## 57       2                   65
    ## 58       3                   62
    ## 59       4                   71
    ## 60       5                   76
    ## 61       0                    0
    ## 62       1                   45
    ## 63       2                   58
    ## 64       3                   64
    ## 65       4                   69
    ## 66       5                   74
    ## 67       0                    0
    ## 68       1                   41
    ## 69       2                   54
    ## 70       3                   63
    ## 71       4                   69
    ## 72       5                   81
    ## 73       0                    0
    ## 74       1                   43
    ## 75       2                   59
    ## 76       3                   58
    ## 77       4                   71
    ## 78       5                   80
    ## 79       0                    0
    ## 80       1                   42
    ## 81       2                   54
    ## 82       3                   65
    ## 83       4                   70
    ## 84       5                   75
    ## 85       0                    0
    ## 86       1                   46
    ## 87       2                   57
    ## 88       3                   66
    ## 89       4                   75
    ## 90       5                   72
    ## 91       0                    0
    ## 92       1                   44
    ## 93       2                   57
    ## 94       3                   70
    ## 95       4                   73
    ## 96       5                   81
    ## 97       0                    0
    ## 98       1                   40
    ## 99       2                   54
    ## 100      3                   55
    ## 101      4                   71
    ## 102      5                   83
    ## 103      0                    0
    ## 104      1                   39
    ## 105      2                   50
    ## 106      3                   70
    ## 107      4                   73
    ## 108      5                   74
    ## 109      0                    0
    ## 110      1                   42
    ## 111      2                   58
    ## 112      3                   70
    ## 113      4                   72
    ## 114      5                   79
    ## 115      0                    0
    ## 116      1                   48
    ## 117      2                   55
    ## 118      3                   72
    ## 119      4                   64
    ## 120      5                   78
    ## 121      0                    0
    ## 122      1                   38
    ## 123      2                   56
    ## 124      3                   62
    ## 125      4                   68
    ## 126      5                   81
    ## 127      0                    0
    ## 128      1                   36
    ## 129      2                   62
    ## 130      3                   66
    ## 131      4                   74
    ## 132      5                   81
    ## 133      0                    0
    ## 134      1                   36
    ## 135      2                   56
    ## 136      3                   72
    ## 137      4                   74
    ## 138      5                   73
    ## 139      0                    0
    ## 140      1                   33
    ## 141      2                   62
    ## 142      3                   67
    ## 143      4                   79
    ## 144      5                   82
    ## 145      0                    0
    ## 146      1                   45
    ## 147      2                   50
    ## 148      3                   66
    ## 149      4                   79
    ## 150      5                   74
    ## 151      0                    0
    ## 152      1                   47
    ## 153      2                   57
    ## 154      3                   73
    ## 155      4                   70
    ## 156      5                   69
    ## 157      0                    0
    ## 158      1                   41
    ## 159      2                   60
    ## 160      3                   69
    ## 161      4                   70
    ## 162      5                   76
    ## 163      0                    0
    ## 164      1                   42
    ## 165      2                   49
    ## 166      3                   73
    ## 167      4                   73
    ## 168      5                   72
    ## 169      0                    0
    ## 170      1                   50
    ## 171      2                   58
    ## 172      3                   65
    ## 173      4                   80
    ## 174      5                   70
    ## 175      0                    0
    ## 176      1                   37
    ## 177      2                   49
    ## 178      3                   68
    ## 179      4                   70
    ## 180      5                   75
    ## 181      0                    0
    ## 182      1                   45
    ## 183      2                   53
    ## 184      3                   70
    ## 185      4                   75
    ## 186      5                   85
    ## 187      0                    0
    ## 188      1                   45
    ## 189      2                   57
    ## 190      3                   63
    ## 191      4                   78
    ## 192      5                   71
    ## 193      0                    0
    ## 194      1                   45
    ## 195      2                   54
    ## 196      3                   70
    ## 197      4                   73
    ## 198      5                   74
    ## 199      0                    0
    ## 200      1                   44
    ## 201      2                   55
    ## 202      3                   64
    ## 203      4                   73
    ## 204      5                   78
    ## 205      0                    0
    ## 206      1                   44
    ## 207      2                   64
    ## 208      3                   74
    ## 209      4                   78
    ## 210      5                   79
    ## 211      0                    0
    ## 212      1                   42
    ## 213      2                   57
    ## 214      3                   70
    ## 215      4                   72
    ## 216      5                   76
    ## 217      0                    0
    ## 218      1                   43
    ## 219      2                   58
    ## 220      3                   67
    ## 221      4                   72
    ## 222      5                   67
    ## 223      0                    0
    ## 224      1                   38
    ## 225      2                   55
    ## 226      3                   67
    ## 227      4                   69
    ## 228      5                   81
    ## 229      0                    0
    ## 230      1                   43
    ## 231      2                   50
    ## 232      3                   69
    ## 233      4                   60
    ## 234      5                   75
    ## 235      0                    0
    ## 236      1                   34
    ## 237      2                   60
    ## 238      3                   61
    ## 239      4                   66
    ## 240      5                   69
    ## 241      0                    0
    ## 242      1                   37
    ## 243      2                   61
    ## 244      3                   71
    ## 245      4                   75
    ## 246      5                   80
    ## 247      0                    0
    ## 248      1                   47
    ## 249      2                   49
    ## 250      3                   63
    ## 251      4                   68
    ## 252      5                   78
    ## 253      0                    0
    ## 254      1                   43
    ## 255      2                   55
    ## 256      3                   55
    ## 257      4                   69
    ## 258      5                   79
    ## 259      0                    0
    ## 260      1                   44
    ## 261      2                   63
    ## 262      3                   62
    ## 263      4                   76
    ## 264      5                   76
    ## 265      0                    0
    ## 266      1                   42
    ## 267      2                   57
    ## 268      3                   70
    ## 269      4                   75
    ## 270      5                   78
    ## 271      0                    0
    ## 272      1                   47
    ## 273      2                   54
    ## 274      3                   68
    ## 275      4                   77
    ## 276      5                   77
    ## 277      0                    0
    ## 278      1                   41
    ## 279      2                   57
    ## 280      3                   69
    ## 281      4                   71
    ## 282      5                   80
    ## 283      0                    0
    ## 284      1                   41
    ## 285      2                   66
    ## 286      3                   68
    ## 287      4                   69
    ## 288      5                   83
    ## 289      0                    0
    ## 290      1                   40
    ## 291      2                   65
    ## 292      3                   62
    ## 293      4                   70
    ## 294      5                   72
    ## 295      0                    0
    ## 296      1                   48
    ## 297      2                   68
    ## 298      3                   70
    ## 299      4                   73
    ## 300      5                   74

``` r
output_frame$fraction_successful <- output_frame$successes_in_hundred / 100
output_frame_summarized <- output_frame %>% group_by(clones) %>% summarize(mean = mean(fraction_successful),
                                                                                       sd = sd(fraction_successful))
```

    ## `summarise()` ungrouping output (override with `.groups` argument)

``` r
Bootstrapped_clone_screening_simulation <- ggplot() + 
  theme_bw() + 
  theme(panel.grid.minor = element_blank(), panel.grid.major.x = element_blank()) +
  labs(x = "Number of clones screened by miniprep + sequencing",
       y = "Fraction of times screening\nthe indicated number of clones\nyielded at least one successful clone\n") +
  scale_y_continuous(limits = c(0,1)) + 
  geom_hline(yintercept = 0) +
  geom_hline(yintercept = 1) +
  geom_line(data = output_frame_summarized, aes(x = clones, y = mean), alpha = 0.4) +
  geom_errorbar(data = output_frame_summarized, aes(x = clones, ymin = mean - sd*1.96, ymax = mean + sd*1.96)
                , width = 0.1, alpha = 0.6) + 
  geom_point(data = output_frame_summarized, aes(x = clones, y = mean))
ggsave(file = "Bootstrapped_clone_screening_simulation.pdf", Bootstrapped_clone_screening_simulation, height = 3, width = 5)
Bootstrapped_clone_screening_simulation
```

![](Cloning_success_rates_files/figure-gfm/Part%20two-2.png)<!-- -->
