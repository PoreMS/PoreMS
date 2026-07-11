# PoreMS Visual Design

## Color Palette

The logo and all visual assets use a strict 3-color blue palette derived from the
[Material Design Blue scale](https://m2.material.io/design/color/the-color-system.html).

| Role | Hex | RGB | Material token |
|---|---|---|---|
| Bond lines | `#90CAF9` | `rgb(144, 202, 249)` | Blue 200 |
| Small/medium nodes, accent text | `#1976D2` | `rgb(25, 118, 210)` | Blue 700 |
| Large anchor nodes, body text | `#0D47A1` | `rgb(13, 71, 161)` | Blue 900 |

### Rationale

- **3 colors maximum** keeps the graphic legible at small sizes (favicon, sidebar logo).
- **Perceptual hierarchy**: bonds are the lightest element (they are connectors, not atoms);
  large anchor nodes are the darkest (they carry the most structural weight).
- **Material Blue scale** provides perceptually uniform steps and is well-known in digital design,
  making the palette easy to extend consistently if needed.
- **Single-hue scheme** avoids colour-meaning confusion — in chemistry, element colours
  (CPK) are standardized. Using a neutral blue avoids implying specific elements.

## Logo Files

| File | Usage |
|---|---|
| `docs/pics/logo.svg` | Square icon (favicon source, app icon) |
| `docs/pics/logo_text.svg` | Horizontal logo with "PoreMS" wordmark |
| `docs/pics/logo_text_sub.svg` | Logo with wordmark + subtitle line |

## Typography

Wordmark uses **Arial / Arial MT** (system sans-serif fallback).
- "**PO**" and "**RE**": Blue 900 (`#0D47A1`)
- "**MS**": Blue 700 (`#1976D2`) — accent color highlights the abbreviation

## Favicon

The favicon is derived from `logo.svg`. At 32 × 32 px the three-color scheme remains
distinguishable because the large anchor node (Blue 900) contrasts clearly against the
bond lines (Blue 200).
