# Copyright (C) 2023  Ferdous,S M <ferdous.csebuet@egmail.com>
# Author: Ferdous,S M <ferdous.csebuet@egmail.com>
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.

import json
import csv
import sys, getopt
from numpy import genfromtxt




injson = open('../data/Ham_N2O.json')
outjson = open('../data/Ham_N2O_part.json','w')

D = json.load(injson)
D_list=list(D.items())

csv_reader = csv.reader(incol,delimiter=' ')

print(csv_reader)

my_data = genfromtxt('../data/Ham_N2O.col',delimiter=' ',dtype="int32")


print(my_data)

partition = {} 

n = len(D_list)

for i in range(n-1):
    X = D_list[i][0]
    c = my_data[i][1]
    c = str(c)
    if c in partition.keys():
        partition[c].append(X)
    else:
        partition[c] = [X]



print(partition)

json.dump(partition,outjson)

