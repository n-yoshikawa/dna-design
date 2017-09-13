#!/usr/bin/env python
# -*- coding: utf-8 -*-

import cgi, cgitb
import dna_design
import constants as c

cgitb.enable()

print('Content-type: text/html')
print('Access-Control-Allow-Origin: *\n')

try:
    form = cgi.FieldStorage()
    text = form.getfirst("text")
    n = form.getfirst("number")
    sequence_list = []
    for amino in text.splitlines():
        sequence_list.append(amino)
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


    for seq, aminos in generated_aminos:
        total = len(aminos)
        inlist = 0
        for amino in aminos:
            if amino in sequence_list:
                inlist += 1
        print('''<table class="table table-striped" style="margin-top:3rem">
          <thead>
            <tr>
              <th>{}から合成されるアミノ酸</th>
              <th>合成対象か？（合成対象の割合: {}/{}）</th>
            </tr>
          </thead>
          <tbody>'''.format(seq, inlist, total))
        for amino in aminos:
            print('<tr>')
            print('<td>{}</td><td>{}</td>'.format(amino, "<span class=\"badge badge-success\">はい</span>" if amino in sequence_list else "<span class=\"badge badge-danger\">いいえ</span>"))
            print('</tr>')
        print('''  </tbody>
        </table>''')
except:
    print('<p>エラーが発生しました。合成するアミノ酸の長さが全て同じかどうか確認してください。</p>')
