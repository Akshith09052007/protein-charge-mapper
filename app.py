from flask import Flask, jsonify, request, render_template
from rna_predictor import predict_structure

app = Flask(__name__)
MAX_WEB_SEQUENCE_LENGTH = 300


@app.get("/")
def home():
    return render_template("index.html")


@app.get("/api/health")
def health():
    return jsonify({
        "status": "ok",
        "service": "RNA Secondary Structure Basic Predictor"
    })


@app.post("/api/predict")
def predict():
    data = request.get_json(silent=True)
    if not isinstance(data, dict):
        return jsonify({"success": False, "error": "Send a JSON object."}), 400

    sequence = data.get("sequence", "")
    if not isinstance(sequence, str):
        return jsonify({"success": False, "error": "Sequence must be text."}), 400
    if len(sequence) > MAX_WEB_SEQUENCE_LENGTH * 2:
        return jsonify({
            "success": False,
            "error": f"Input is too long. Maximum sequence length is {MAX_WEB_SEQUENCE_LENGTH} nt."
        }), 400

    try:
        min_loop = data.get("min_loop", 3)
        allow_wobble = data.get("allow_wobble", True)
        result = predict_structure(sequence, min_loop, allow_wobble)
        if result["length"] > MAX_WEB_SEQUENCE_LENGTH:
            return jsonify({
                "success": False,
                "error": f"Sequence is too long. Maximum sequence length is {MAX_WEB_SEQUENCE_LENGTH} nt."
            }), 400
        return jsonify({"success": True, "result": result})
    except (ValueError, TypeError) as error:
        return jsonify({"success": False, "error": str(error)}), 400


if __name__ == "__main__":
    app.run(host="127.0.0.1", port=5000, debug=False)
