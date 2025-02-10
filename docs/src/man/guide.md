# Guide

You can create a simple XML by running
```julia
using PhysiCellECMCreator
path_to_folder = "template_folder"
createICECMXMLTemplate(path_to_folder) # will create ecm.xml inside the folder "template_folder"
```

Then, you can edit the XML file as needed.

## Layers
The ECM is divided into layers.
The layers decide precedence of ECM values in the voxels.
The lower layers (lower `ID` values) are drawn first, and the higher layers are drawn on top, overwriting the lower layers.
Within each layer, patches cannot overlap.

## Patches
Each layer can have multiple patches.
The patches can be of different types, selected from `everywhere`, `ellipse`, and `elliptical_disc`.
See [Example XML](@ref) for an example of each type.

### `everywhere`
The `everywhere` type is used to set the ECM values for the entire domain.
This is useful for the base layer, i.e., `ID="1"`.
You may set the `density`, `orientation`, and `anisotropy` for the entire domain.
The only `orientation` option is `random`.
A sample `everywhere` patch is shown below.
```xml
<patch ID="1" type="everywhere">
    <density>0.4</density>
    <orientation>random</orientation>
    <anisotropy>0.3</anisotropy>
</patch>
```

### `ellipse`
The `ellipse` type is used to set the ECM values on an ellipse.
The `x0` and `y0` elements set the center of the ellipse.
The `a` and `b` elements set the axes in the x and y directions, respectively.
Thus, the equation of the ellipse considering only these parameters is

```math
\begin{aligned}
\frac{(x-x_0)^2}{a^2} + \frac{(y-y_0)^2}{b^2} = 1
\end{aligned}
```

The `thickness` element can be set to a positive value to create a ring of ECM that begins at the ellipse and extends outwards to the ellipse defined by

```math
\begin{aligned}
\frac{(x-x_0)^2}{(a+thickness)^2} + \frac{(y-y_0)^2}{(b+thickness)^2} = 1
\end{aligned}
```

The `rotation` element sets the counter-clockwise rotation of the ellipse.
By default, `rotation` is assumed to be in degrees.
If the content of the element includes "pi", then the units are assumed to be in radians.
Furthermore, the units can be explicitly set with the `units` attribute (see below).
So long as the attribute satisfies `startswith(units, "rad")`, the units will be assumed to be in radians.

The `orientation` element can be set to `parallel`, `perpendicular`, or `random`.
For `parallel`, the orientation will follow the shape of the ellipse.
For `perpendicular`, the orientation will be perpendicular to the shape of the ellipse.
A sample `ellipse` patch is shown below.
```xml
<patch ID="1" type="ellipse">
    <density>0.5</density>
    <orientation>parallel</orientation>
    <anisotropy>1.0</anisotropy>
    <x0>150.0</x0>
    <y0>200.0</y0>
    <a>100.0</a>
    <b>80.0</b>
    <rotation units="rad">2pi/3</rotation>
    <thickness>100.0</thickness>
</patch>
```

### `elliptical_disc`
The `elliptical_disc` type is used to set the ECM values on an elliptical disc.
It is often used in tandem with an `ellipse` patch to fill in the center of the ellipse.
The same parameters as the `ellipse` patch are used, except the `thickness` parameter is not used.
The only `orientation` option is `random`.
A sample `elliptical_disc` patch is shown below.
```xml
<patch ID="2" type="elliptical_disc">
    <density>0.3</density>
    <orientation>random</orientation>
    <anisotropy>0.3</anisotropy>
    <x0>150.0</x0>
    <y0>200.0</y0>
    <a>100.0</a>
    <b>80.0</b>
    <rotation units="rad">2pi/3</rotation>
</patch>
```

### `ellipse_with_shell`
The `ellipse_with_shell` type is used to set the ECM values on an ellipse with a shell, i.e. a combination of `ellipse` and `elliptical_disc`.
The location of the ellipse is set exactly as the `ellipse` and `elliptical_disc` patches.
However, the `density`, `orientation`, and `anisotropy` values are set for the `interior` and `shell` in separate child elements with those tags.
The same orientation restrictions for `ellipse` apply to the `shell`.
Similarly, `elliptical_disc` orientation restrictions apply to the `interior`.
A sample `ellipse_with_shell` patch is shown below.
```xml
<patch ID="3" type="ellipse_with_shell">
    <x0>150.0</x0>
    <y0>200.0</y0>
    <a>100.0</a>
    <b>80.0</b>
    <rotation units="rad">2pi/3</rotation>
    <thickness units="px">100.0</thickness>
    <interior>
        <density>0.5</density>
        <orientation>random</orientation>
        <anisotropy>1.0</anisotropy>
    </interior>
    <shell>
        <density>0.3</density>
        <orientation>parallel</orientation>
        <anisotropy>0.3</anisotropy>
    </shell>
</patch>
```
An additional `exterior` element can be added analagously to the `interior` and `shell` elements.
This will define the ECM values outside the shell.
The only `orientation` option is `random`.
Caution: using the `exterior` option means the `ellipse_with_shell` patch is setting ECM values throughout the entire domain.
Any previous layer will be overwritten.
Any patches in the same layer will necessarily overlap, causing an error.

## Example XML

A full example XML is shown below.
```xml
<?xml version="1.0"?>
<ic_ecm>
    <layer ID="1">
        <patch_collection type="everywhere">
            <patch ID="1">
                <density>0.4</density>
                <orientation>random</orientation>
                <anisotropy>0.3</anisotropy>
            </patch>
        </patch_collection>
        ...
    </layer>
    ...
    <layer ID="3">
        <patch_collection type="ellipse">
            <patch ID="1">
                <density>0.5</density>
                <orientation>parallel</orientation>
                <anisotropy>1.0</anisotropy>
                <x0>150.0</x0>
                <y0>200.0</y0>
                <a>100.0</a>
                <b>80.0</b>
                <rotation units="rad">2pi/3</rotation>
                <thickness units="px">100.0</thickness>
            </patch>
        </patch_collection>
        <patch_collection type="elliptical_disc">
            <patch ID="1">
                <density>0.3</density>
                <orientation>random</orientation>
                <anisotropy>0.3</anisotropy>
                <x0>150.0</x0>
                <y0>200.0</y0>
                <a>100.0</a>
                <b>80.0</b>
                <rotation units="rad">2pi/3</rotation>
            </patch>
        </patch_collection>
        <patch_collection type="ellipse_with_shell">
            <patch ID="1">
                <x0>150.0</x0>
                <y0>200.0</y0>
                <a>100.0</a>
                <b>80.0</b>
                <rotation units="rad">2pi/3</rotation>
                <thickness units="px">100.0</thickness>
                <interior>
                    <density>0.5</density>
                    <orientation>random</orientation>
                    <anisotropy>1.0</anisotropy>
                </interior>
                <shell>
                    <density>0.3</density>
                    <orientation>parallel</orientation>
                    <anisotropy>0.3</anisotropy>
                </shell>
            </patch>
        </patch_collection>
        ...
    </layer>
</ic_ecm>
```