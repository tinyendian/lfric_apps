##############################################################################
# (c) Crown copyright 2024 Met Office. All rights reserved.
# The file LICENCE, distributed with this code, contains details of the terms
# under which the code may be used.
##############################################################################
# Wrapper script to build the required PSyAD targets.
# The PSyAD kernels are built using three stages.
# 1) Pre-patch: copies (and patches) tangent linear kernels from their base dir to PSYAD_WDIR.
# 2) PSyAD: generates adjoint kernels and adjoint test algorithms from the pre-patch stage tangent linear kernels.
# 3) Post-patch: copies (and patches) adjoint kernels and adjoint test algorithms from PSYAD_WDIR to WORKING_DIR.

# TODO: Setup the psyad command to use '--config=/path/to/psyclone.cfg'.
# This requires the PSyclone issue https://github.com/stfc/PSyclone/issues/2826 to be fixed.

# List of targets needed by this script.
include $(ADJOINT_BUILD)/psyad_targets.mk

all: $(POST_PATCH_TARGETS) $(PSYAD_TARGETS) $(PRE_PATCH_TARGETS)

# List of rules to build the adjoint kernels and adjoint test algorithms.
include $(ADJOINT_BUILD)/pre_patch.mk
include $(ADJOINT_BUILD)/psyad.mk
include $(ADJOINT_BUILD)/post_patch.mk

$(DIRECTORIES):
	mkdir -p $@
