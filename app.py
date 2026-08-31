from flask import Flask, jsonify, request, render_template
from rna_predictor import predict_structure

app = Flask(__name__)

@app.get("/")
def home():
    return render_template("index.html")

@app.get("/api/health")
def health():
    return jsonify({"status": "ok", "service": "RNA Secondary Structure Predictor"})

@app.post("/api/predict")
def predict():
    data = request.get_json(silent=True) or {}
    sequence = data.get("sequence", "")
    try:
        result = predict_structure(
            sequence,
            int(data.get("min_loop", 3)),
            bool(data.get("allow_wobble", True)),
        )
        return jsonify({"success": True, "result": result})
    except (ValueError, TypeError) as exc:
        return jsonify({"success": False, "error": str(exc)}), 400

if __name__ == "__main__":
    app.run(host="0.0.0.0", port=5000, debug=True)
