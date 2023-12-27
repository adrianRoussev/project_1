from flask import Flask, jsonify, request
from flask_restful import Api, Resource
from convert_ord_reaction_data import ReactionSimilarity  

app = Flask(__name__)
api = Api(app)

class ReactionSimilarityResource(Resource):
    def post(self):
        data = request.get_json()
        file = data.get('file')
        input_smiles = data.get('input_smiles')
        cutoff_score = data.get('cutoff_score')
        similarity_type = data.get('similarity_type')
        params = data.get('params')
        function_name = data.get('function_name')  # Add a new field for function_name

        reaction_similarity = ReactionSimilarity(file, input_smiles, cutoff_score, similarity_type, params)

        # Use getattr to dynamically call the method based on function_name
        try:
            method_to_call = getattr(reaction_similarity, function_name)
        except AttributeError:
            return jsonify({'error': 'Invalid function_name'})

        # Call the method and return the result
        result = method_to_call()
        return jsonify(result)

api.add_resource(ReactionSimilarityResource, '/reaction-similarity')

if __name__ == '__main__':
    app.run(debug=True, port=9000)

    # curl -X POST -H "Content-Type: application/json" -d '{"file": "/path/to/your/file","input_smiles": "[CH3:1][NH:2][CH3:3].Cl[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1>>[CH3:1][N:2]([CH3:3])[CH2:5][CH2:6][CH:7]1[CH2:13][N:12]([CH3:14])[C:11](=[S:15])[C:10]2[CH:16]=[C:17]([N+:20]([O-:22])=[O:21])[CH:18]=[CH:19][C:9]=2[O:8]1","cutoff_score": 0.5,"similarity_type": "total similarity","params": ["identifiers", "reaction_id"],"function_name": "calculate_similarity_from_file"}' http://localhost:5000/reaction-similarity