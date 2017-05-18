import csv

with open('infile','r'), open ('outfile','w') as fin, fout:
    writer = csv.writer(fout, delimiter=' ')            
    for row in csv.reader(fin, delimiter=' '):
        if row[2] == 'Central':
             writer.writerow(row)