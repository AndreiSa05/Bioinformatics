
import os
import io
from contextlib import redirect_stdout

from dotenv import load_dotenv
from openai import OpenAI

load_dotenv()

with redirect_stdout(io.StringIO()):
    import ex1

virus_name = ex1.virus_name
top10 = ex1.top10

client = OpenAI(api_key=os.getenv("OPENAI_API_KEY"))

prompt = f"""
Virus identified from FASTA:
{virus_name}

Top 10 nucleotide triplets:
{top10}

Provide only practical food and beverage recommendations
for someone recovering from an infection with this virus.

Output exactly two sections:

FOODS:
- List recommended foods.

DRINKS:
- List recommended drinks and general hydration advice.

Rules:
- No introduction or conclusion.
- No genome or codon analysis in the output.
- No explanations unless absolutely necessary.
- Do not derive food recommendations from codon frequencies.
- Recommend foods based on established medical knowledge.
- Do not invent treatments or claim food cures a virus.
- If this virus is not known to infect humans, or there is
  no applicable evidence, do not invent a recovery diet.
  Instead, state briefly that no virus-specific
  recommendations are available.
"""

response = client.responses.create(
    model="gpt-4.1-mini",
    input=prompt
)

print(response.output_text)