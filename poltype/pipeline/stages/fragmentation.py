"""
poltype.pipeline.stages.fragmentation – molecular fragmentation stage.

Decomposes the current molecule into fragments centred on rotatable
bonds that require torsion scanning.  The :class:`FragmentResult` is
stored as the ``"fragment_result"`` artifact for consumption by
downstream stages (primarily torsion fitting).

This stage is skipped when ``config.use_fragmentation`` is ``False``,
``config.database_match_only`` is ``True``, or the pipeline is
running in dry-run mode.
"""

from __future__ import annotations

import logging
from typing import TYPE_CHECKING, Optional

from poltype.fragmentation.fragment import FragmentResult
from poltype.fragmentation.fragmenter import (
    Fragmenter,
    RuleBasedFragmenter,
    WBOFragmenter,
)
from poltype.pipeline.stage import Stage, StageResult, StageStatus

if TYPE_CHECKING:
    from poltype.pipeline.context import PipelineContext

logger = logging.getLogger(__name__)


class FragmentationStage(Stage):
    """Break a molecule into fragments for efficient torsion scanning.

    Parameters
    ----------
    fragmenter:
        A :class:`Fragmenter` instance.  If ``None``, a
        :class:`RuleBasedFragmenter` is built during execution.
    """

    def __init__(self, fragmenter: Optional[Fragmenter] = None) -> None:
        super().__init__(name="fragmentation")
        self._fragmenter = fragmenter

    def execute(self, context: "PipelineContext") -> StageResult:
        """Fragment the molecule around rotatable bonds.

        Uses the ``db_match_result`` artifact (if available) to
        determine which rotatable bonds still need torsion scanning,
        and falls back to the full list of rotatable bonds when no
        database match result is present.

        Returns
        -------
        StageResult
            ``COMPLETED`` with a ``fragment_result`` artifact on
            success, ``FAILED`` if there is no molecule.
        """
        mol = context.molecule
        if mol is None:
            return StageResult(
                status=StageStatus.FAILED,
                message="No molecule in context",
            )

        # Determine which rotatable bonds need scanning.
        rotatable = mol.rotatable_bonds
        db_result = context.get_artifact("db_match_result")
        if db_result is not None:
            # Filter to bonds whose atoms are not fully matched.
            matched = set(db_result.matched_torsions)
            needs_scan = [
                bond for bond in rotatable
                if not self._torsion_covered(bond, matched)
            ]
        else:
            needs_scan = list(rotatable)

        if not needs_scan:
            logger.info("No rotatable bonds need scanning; skipping fragmentation")
            result = FragmentResult(
                parent_molecule_name=mol.name,
                fragmentation_method="none",
            )
            context.set_artifact("fragment_result", result)
            return StageResult(
                status=StageStatus.COMPLETED,
                message="No rotatable bonds to fragment",
                artifacts={"fragment_result": result},
            )

        # Build fragmenter lazily from config if not injected.
        fragmenter = self._fragmenter
        if fragmenter is None:
            fragmenter = RuleBasedFragmenter()

        frag_result: FragmentResult = fragmenter.fragment(
            mol,
            needs_scan,
        )
        context.set_artifact("fragment_result", frag_result)

        logger.info(
            "Fragmentation complete: %d fragment(s) for %d bond(s)",
            frag_result.num_fragments,
            len(needs_scan),
        )
        return StageResult(
            status=StageStatus.COMPLETED,
            message=(
                f"{frag_result.num_fragments} fragment(s) "
                f"for {len(needs_scan)} rotatable bond(s)"
            ),
            artifacts={"fragment_result": frag_result},
        )

    def should_skip(self, context: "PipelineContext") -> bool:
        """Skip when fragmentation is disabled, database-only, or dry-run."""
        return (
            not context.config.use_fragmentation
            or context.config.database_match_only
            or context.config.dry_run
        )

    @staticmethod
    def _torsion_covered(
        bond: tuple,
        matched_torsions: set,
    ) -> bool:
        """Return ``True`` if the rotatable *bond* is covered by a matched torsion.

        A bond ``(i, j)`` is considered covered if *any* matched
        torsion tuple contains both *i* and *j* as the central pair
        (positions 1 and 2 in the 4-tuple).
        """
        for torsion in matched_torsions:
            if len(torsion) >= 4:
                central = (torsion[1], torsion[2])
                if (bond[0], bond[1]) == central or (bond[1], bond[0]) == central:
                    return True
        return False
