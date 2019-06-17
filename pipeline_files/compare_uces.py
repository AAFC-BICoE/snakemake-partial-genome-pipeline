"""
Author: Jackson Eyres
Copyright: Government of Canada
License: MIT
"""

import os
import argparse
import csv
import sqlite3
import operator
from itertools import groupby


def main():
    parser = argparse.ArgumentParser(description='Compares UCEs found in Phyluce generated sqlite database')
    parser.add_argument('-db', nargs="+",
                        help='SQLite Database files', required=True)

    args = parser.parse_args()
    #print(args)
    csv_files = convert_sqlite_csv(args.db)
    analyze_csv(csv_files)

def convert_sqlite_csv(database_list):
    """
    Takes a Phyluce generated probe.matches.sqlite database, and generates a match_map table in CSV format for easier parsing.
    :param database_directory:
    :return:
    """
    csv_files = []
    for database in database_list:
        database_directory = os.path.split(database)[0]
        file_name = os.path.split(database)[1]
        with sqlite3.connect(database) as connection:
            csv_file = os.path.join(database_directory, file_name.replace("sqlite", "csv"))
            csvWriter = csv.writer(open(csv_file, "w"))
            c = connection.cursor()
            for tablename in c.execute('SELECT name FROM sqlite_master WHERE type="table";'):
                if tablename[0].startswith('match_map'):
                    tablename = str(tablename[0])
                    query = "Select * FROM " + tablename
                    query = str(query)
                    data = c.execute(query)
                    # print(c.description)
                    names = list(map(lambda x: x[0], c.description))
                    # print(names)
                    csvWriter.writerow(names)
                    csvWriter.writerows(data)
        csv_files.append(csv_file)
    return csv_files


def analyze_csv(csv_files):
    uce_dict = {}

    for csv_file in csv_files:
        with open (csv_file) as f:
            lines = f.readlines()
            lines.pop(0)
            print("Number of Samples in {}: {}".format(csv_file, len(lines[0].split(","))-1))
            for line in lines:
                split = line.rstrip().split(",")
                count = -1
                for item in split:
                    if item is not "":
                        count += 1

                if split[0] in uce_dict:
                    uce_dict[split[0]] += count
                else:
                    uce_dict[split[0]] = count

    sorted_x = sorted(uce_dict.items(), key=operator.itemgetter(1))
    for tuple in sorted_x:
        print(tuple[0],tuple[1])

    X2 = [list(group) for key, group in groupby(sorted_x, operator.itemgetter(1))]
    print(X2)
    for item in X2:
        print(item[0][1], len(item))
if __name__ == "__main__":
    main()
