## Created: April 18, 2024
## Updated: April 19, 2024
## Author(s): Todd Lenz, tlenz001@ucr.edu

# Define variables
INSTALL_SCRIPT = INSTALL
ENVIRONMENT_FILE = environment.yml

# Default target
all: install

# Install package
install:
	source ./$(INSTALL_SCRIPT) 

# Clean package
clean:
	rm -rf $(ENVIRONMENT_FILE).lock
	rm -rf $(ENVIRONMENT_FILE).bak
	conda env remove -y -n SimpleSeq

# Uninstall package
uninstall: clean
	echo "Do you want to uninstall SimpleSeq? (y/n)"
	read -r ans
	if [[ "$ans" = "y" ]]; then
		sed -i '/SimpleSeq/d' ~/.bashrc
		source ~/.bashrc
	else
		echo "Uninstallation aborted."
	fi
