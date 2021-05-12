install:
	pip3 -v install -v .

uninstall:
	pip3 -v uninstall -y manifold

test:
	python3 -m unittest discover -v tests

editable:
	pip3 -v install -e . && cp build/*/compile_commands.json build/
