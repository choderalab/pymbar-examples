import numpy
import re
import os

#===================================================================================================
# FUNCTIONS: The unix-like helpers.
#===================================================================================================

def wcPy(f):
   """Count up lines in file 'f'."""
   if not type(f) is file:
      with open(f, 'r') as f:
         return wcPy(f)
   return sum(1 for l in f)

def trPy(s, l='[,\\\\"/()-]', char=' '):
   """In string 's' replace all the charachters from 'l' with 'char'."""
   return re.sub(l, char, s)

def grepPy(f, s):
   """From file 'f' extract the (first occurence of) line that contains string 's'."""
   if not type(f) is file:
      with open(f, 'r') as f:
         return grepPy(f, s)
   for line in f:
      if s in line:
         return line
   return ''

def tailPy(f, nlines, lenb=1024):
   if not type(f) is file:
      with open(f, 'r') as f:
         return tailPy(f, nlines, lenb)
   f.seek(0, 2)
   sizeb = f.tell()
   n_togo = nlines
   i = 1
   excerpt = []
   while n_togo > 0 and sizeb > 0:
      if (sizeb - lenb > 0):
         f.seek(-i*lenb, 2)
         excerpt.append(f.read(lenb))
      else:
         f.seek(0,0)
         excerpt.append(f.read(sizeb))
      ll = excerpt[-1].count('\n')
      n_togo -= ll
      sizeb -= lenb
      i += 1
   return ''.join(excerpt).splitlines()[-nlines:]


