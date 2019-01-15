import sys

from flask import Flask, request, render_template, jsonify

# import dna_design
import clustering

app = Flask(__name__)


@app.route('/')
def index():
    return render_template('index.html')


@app.route("/api", methods=['POST'])
def hello():
    sequence_list = []
    for amino in request.json["sequences"].splitlines():
        sequence_list.append(amino)
    n = -1 if request.json["type"] == "auto" else int(request.json["number"])
    scale = float(request.json["scale"])
    # result_base, generated_aminos = dna_design.optimize(sequence_list, n)
    result_base, generated_aminos = clustering.design(sequence_list, n, scale)

    results = []
    # print("sequence_list:", sequence_list)
    size = 0
    for base, amino in zip(result_base, generated_aminos):
        # print(base)
        # print(amino)
        results.append({
            "base": base,
            "amino": [{"codon": a[0], "seq": a[1], "isTarget": a[1] in sequence_list} for a in amino],
            "count": sum([a[1] in sequence_list for a in list(set(amino))])})
        size += len(amino)
        print(len(amino))
    data = {
        "results": results,
        "sequences": sequence_list,
        "size": size
    }
    return jsonify(data)


class Unbuffered(object):
    def __init__(self, stream):
        self.stream = stream

    def write(self, data):
        self.stream.write(data)
        self.stream.flush()

    def writelines(self, datas):
        self.stream.writelines(datas)
        self.stream.flush()

    def __getattr__(self, attr):
        return getattr(self.stream, attr)


sys.stdout = Unbuffered(sys.stdout)

if __name__ == '__main__':
    app.run(debug=True, port=5555)
