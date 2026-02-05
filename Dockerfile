# 1️ Use official Python image
FROM python:3.10-slim

# ️⃣ Set working directory inside container
WORKDIR /app

# 3️ Copy project files into container
COPY . /app

# 4 Install Python dependencies
RUN pip install --no-cache-dir -r requirements.txt

# 5️ Command to run your code
CMD ["python", "gene2.py"]
