#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cgi, cgitb
import dna_design
import constants as c



print "Optimal sequence"

cgitb.enable()

form = cgi.FieldStorage()
text = form.getfirst("text")
n = form.getfirst("number")
sequence_list = []
for amino in text.splitlines():
    sequence_list.append(amino)
print('Content-type: text/html\nAccess-Control-Allow-Origin: *\n')

try:
    result_base, generated_aminos = dna_design.optimize(sequence_list, int(n))
    print('''<table class="table table-striped">
      <thead>
        <tr>
          <th>設計した塩基配列</th>
        </tr>
      </thead>
      <tbody>''')
    for b in result_base:
        print('<tr>')
        print('<td>{}</td>'.format(b))
        print('</tr>')
    print('''  </tbody>
    </table>''')

    print('''<table class="table table-striped">
      <thead>
        <tr>
          <th>合成されるアミノ酸</th>
          <td>合成対象か？</td>
        </tr>
      </thead>
      <tbody>''')

    for g in generated_aminos:
        print('<tr>')
        print('<td>{}</td><td>{}</td>'.format(g, "<span class=\"badge badge-success\">はい</span>" if g in sequence_list else "<span class=\"badge badge-danger\">いいえ</span>"))
        print('</tr>')
    print('''  </tbody>
    </table>''')
except:
    print('<p>エラーが発生しました。合成するアミノ酸の長さが全て同じかどうか確認してください。</p>')
