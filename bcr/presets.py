"""Storage for assembly layout presets.

Built-in presets live in `assemble.PRESETS` and cannot be edited or removed.
Presets the user builds are kept in `data/presets.json` and merged over the
built-ins by key.

JSON rather than CSV: a layout is an ordered list of slot objects, and nesting
that into a CSV cell would mean inventing an encoding for it. The parts and
antibody tables stay CSV because their rows are flat.
"""

import json
import os
import tempfile
from pathlib import Path

from . import assemble

DEFAULT_DATA_DIR = Path(__file__).resolve().parents[1] / "data"

# Keys that always come from code, so a saved preset can never shadow one.
RESERVED_KEYS = frozenset(assemble.PRESETS)


def slugify(label: str) -> str:
    """A key from a label: lowercase, words joined by underscores."""
    cleaned = [c.lower() if (c.isalnum() or c in " -_") else " " for c in label]
    parts = "".join(cleaned).replace("-", " ").replace("_", " ").split()
    return "_".join(parts) or "preset"


class PresetStore:
    def __init__(self, data_dir=None):
        self.path = Path(data_dir or DEFAULT_DATA_DIR) / "presets.json"

    # --- reading ---------------------------------------------------------

    def read_user(self) -> dict:
        """User presets as stored. A damaged file yields nothing rather than
        taking the whole app down; the built-ins still work."""
        if not self.path.exists():
            return {}
        try:
            with open(self.path, encoding="utf-8") as fh:
                data = json.load(fh)
        except (OSError, json.JSONDecodeError):
            return {}
        if not isinstance(data, dict):
            return {}
        return {
            key: value
            for key, value in data.items()
            if isinstance(value, dict) and isinstance(value.get("slots"), list)
        }

    def all(self) -> dict:
        """Built-ins plus user presets, each tagged with `builtin`."""
        merged = {}
        for key, value in assemble.PRESETS.items():
            merged[key] = {
                "label": value["label"],
                "description": value.get("description", ""),
                "slots": value["slots"],
                "builtin": True,
                # Express layouts are the ones offered on the back of the
                # assembly card: they place whole chains, so they need nothing
                # from the user but the signal peptides and linkers.
                "express": bool(value.get("express")),
            }
        for key, value in self.read_user().items():
            if key in RESERVED_KEYS:
                continue  # a built-in key is never shadowed
            merged[key] = {
                "label": value.get("label") or key,
                "description": value.get("description", ""),
                "slots": value["slots"],
                "builtin": False,
                "express": False,
                "updated_utc": value.get("updated_utc", ""),
            }
        return merged

    def slots_for(self, key: str):
        """The slot list for a preset key, or None if there is no such key."""
        entry = self.all().get(key)
        return entry["slots"] if entry else None

    # --- writing ---------------------------------------------------------

    def _write(self, data: dict) -> None:
        self.path.parent.mkdir(parents=True, exist_ok=True)
        fd, tmp = tempfile.mkstemp(
            dir=self.path.parent, prefix=self.path.name, suffix=".tmp"
        )
        try:
            with os.fdopen(fd, "w", encoding="utf-8") as fh:
                json.dump(data, fh, indent=2, sort_keys=True)
                fh.write("\n")
            os.replace(tmp, self.path)
        except BaseException:
            if os.path.exists(tmp):
                os.unlink(tmp)
            raise

    def save(self, key, label, description, slots, updated_utc) -> str:
        """Insert or replace a user preset. Returns "added" or "updated".

        The caller validates `slots` with assemble.validate_slots first.
        """
        if key in RESERVED_KEYS:
            raise ValueError(
                f"'{key}' is a built-in preset and cannot be overwritten; "
                f"give this one a different name"
            )
        data = self.read_user()
        outcome = "updated" if key in data else "added"
        data[key] = {
            "label": label,
            "description": description,
            "slots": slots,
            "updated_utc": updated_utc,
        }
        self._write(data)
        return outcome

    def rename_part(self, old: str, new: str) -> list[str]:
        """Point every user preset slot naming `old` at `new` instead.

        Renaming a part would otherwise silently break the layouts built on
        it: a slot keeps the old name, and the preset becomes unbuildable with
        a "missing part" alert. Returns the keys that were rewritten.
        """
        data = self.read_user()
        touched = []
        for key, value in data.items():
            slots = value.get("slots") or []
            if not any(
                isinstance(s, dict) and s.get("part") == old for s in slots
            ):
                continue
            for slot in slots:
                if isinstance(slot, dict) and slot.get("part") == old:
                    slot["part"] = new
            touched.append(key)
        if touched:
            self._write(data)
        return sorted(touched)

    def delete(self, key) -> bool:
        if key in RESERVED_KEYS:
            raise ValueError(f"'{key}' is a built-in preset and cannot be removed")
        data = self.read_user()
        if key not in data:
            return False
        del data[key]
        self._write(data)
        return True
