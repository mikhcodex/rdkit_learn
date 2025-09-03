from __future__ import annotations

from typing import Optional, Sequence, Tuple, Union

from rdkit import Chem
from rdkit.Chem import Draw, rdDepictor


MolLike = Union[Chem.Mol, str]


def show_molecule(
    mol: MolLike,
    title: Optional[str] = None,
    size: Tuple[int, int] = (400, 400),
    kekulize: bool = True,
    highlight_atoms: Optional[Sequence[int]] = None,
    highlight_bonds: Optional[Sequence[int]] = None,
    show_atom_indices: bool = False,
    show_bond_indices: bool = False,
    block: bool = True,
) -> None:
    """
    Display an RDKit molecule in a separate window, similar to matplotlib plots.

    - Accepts an RDKit Mol or a SMILES string.
    - Uses matplotlib for display if available; falls back to PIL image viewer otherwise.

    Parameters
    - mol: RDKit Mol or SMILES string.
    - title: Optional title for the window/figure.
    - size: (width, height) in pixels for the rendered image.
    - kekulize: Whether to kekulize the depiction when possible.
    - highlight_atoms: Optional list of atom indices to highlight.
    - highlight_bonds: Optional list of bond indices to highlight.
    - show_atom_indices: If True, draw atom indices next to atoms.
    - show_bond_indices: If True, draw bond indices next to bonds.
    - block: If True, blocks on `plt.show()`; if False, shows non-blocking.

    Notes
    - matplotlib is an optional dependency. If not installed, a PIL-based viewer is used.
    - 2D coordinates are computed if missing.
    """

    # Convert SMILES to Mol if needed
    if isinstance(mol, str):
        rd_mol = Chem.MolFromSmiles(mol)
        if rd_mol is None:
            raise ValueError("Invalid SMILES string provided.")
    else:
        rd_mol = mol
        if rd_mol is None:
            raise ValueError("None provided as molecule.")

    # Compute 2D coordinates for depiction if not present
    try:
        rdDepictor.Compute2DCoords(rd_mol)
    except Exception:
        # If this fails, RDKit will still try to depict with whatever is available
        pass

    # Render to a PIL image
    if show_atom_indices or show_bond_indices:
        # Use rdMolDraw2D to include indices
        import io
        from rdkit.Chem.Draw import rdMolDraw2D
        from PIL import Image as PILImage  # type: ignore

        drawer = rdMolDraw2D.MolDraw2DCairo(int(size[0]), int(size[1]))
        opts = drawer.drawOptions()
        opts.addAtomIndices = bool(show_atom_indices)
        opts.addBondIndices = bool(show_bond_indices)

        rdMolDraw2D.PrepareAndDrawMolecule(
            drawer,
            rd_mol,
            highlightAtoms=list(highlight_atoms) if highlight_atoms else None,
            highlightBonds=list(highlight_bonds) if highlight_bonds else None,
            kekulize=bool(kekulize),
        )
        drawer.FinishDrawing()
        png = drawer.GetDrawingText()
        img = PILImage.open(io.BytesIO(png))
    else:
        img = Draw.MolToImage(
            rd_mol,
            size=size,
            kekulize=kekulize,
            highlightAtoms=list(highlight_atoms) if highlight_atoms else None,
            highlightBonds=list(highlight_bonds) if highlight_bonds else None,
        )

    # Try matplotlib first for a familiar plotting window
    try:
        import matplotlib.pyplot as plt  # type: ignore

        # Map pixel size to figure size at 100 dpi
        fig_w, fig_h = size[0] / 100.0, size[1] / 100.0
        plt.figure(figsize=(fig_w, fig_h), dpi=100)
        plt.imshow(img)
        plt.axis("off")
        if title:
            plt.title(title)
        if block:
            plt.show()
        else:
            plt.show(block=False)
            # Brief pause to ensure window appears in non-blocking mode
            plt.pause(0.001)
        return
    except ImportError:
        # Fall back to PIL's default image viewer
        if title:
            try:
                # Some environments support setting window title via PIL show
                img.info["title"] = title
            except Exception:
                pass
        img.show()
