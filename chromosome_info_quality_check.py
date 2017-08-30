import pandas
import matplotlib.pyplot as plt

filenameList = ['PTBP1_r1_chr1_test.bed', 'PTBP1_r2_chr1_test.bed']
groupcount = {}
main_groupcount = []
in_both = [0 for i in range(5)]
not_in_both = [0 for i in range(5)]

for filename in filenameList:
    df = pandas.read_table(filename, sep='\t', names=['chrom', 'start', 'stop', 'name', 'score', 'strand'])
    grouped = df.groupby(['chrom', 'start', 'strand'])
    groupcount = grouped.size().to_dict()
    main_groupcount.append(groupcount)

file_1 = set(main_groupcount[0].keys())
file_2 = set(main_groupcount[1].keys())

difference_1 = file_1.difference(file_2)
difference_2 = file_2.difference(file_1)
common_both = set(main_groupcount[0].keys()).intersection(main_groupcount[1].keys())

Dict_1 = main_groupcount[0]
Dict_2 = main_groupcount[1]

for x in common_both:
    count = Dict_1[x] + Dict_2[x]
    in_both[min(count - 1, 4)] += 1

for y in difference_1:
    count_1 = Dict_1[y]
    not_in_both[min(count_1 - 1, 4)] += 1

for z in difference_2:
    count_2 = Dict_2[z]
    not_in_both[min(count_2 - 1, 4)] += 1

print("In both files: ", in_both)
print("In either file :", not_in_both)

sum_both = ([x + y for x, y in zip(in_both, not_in_both)])
Ratio_both = ([x / y for x, y in zip(in_both, sum_both)])

col_label = ['container_1', 'container_2', 'container_3', 'container_4', 'container_5+']
chart = pandas.DataFrame(dict(zip(col_label, [Ratio_both])))
plt.plot(chart)
plt.show()
