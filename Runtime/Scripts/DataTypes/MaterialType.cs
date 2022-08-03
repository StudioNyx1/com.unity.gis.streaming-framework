﻿using System;

namespace Unity.Geospatial.Streaming
{
    /// <summary>
    /// Determine which <see cref="UGMaterial"/> the <see cref="UGMaterialFactory"/> should instantiate.
    /// </summary>
    public readonly struct MaterialType : IEquatable<MaterialType>
    {
        /// <summary>
        /// Determine if the material is lit or unlit
        /// </summary>
        public readonly MaterialLighting Lighting;

        /// <summary>
        /// Determine if the material is opaque, alpha tested, or transparent
        /// </summary>
        public readonly MaterialAlphaMode AlphaMode;

        /// <summary>
        /// Default constructor.
        /// </summary>
        /// <param name="lighting">
        /// Specify whether or not the <see href="https://docs.unity3d.com/ScriptReference/Material.html">material</see>
        /// will be affected by the scene lighting or not.
        /// </param>
        /// <param name="alphaMode">
        /// Specify how decode <see href="https://docs.unity3d.com/ScriptReference/Material.html">Materials</see>
        /// with transparency directives.
        /// </param>
        public MaterialType(MaterialLighting lighting, MaterialAlphaMode alphaMode)
        {
            Lighting = lighting;
            AlphaMode = alphaMode;
        }

        /// <summary>
        /// Get if this instance is the same as the given <paramref name="obj"/>.
        /// </summary>
        /// <param name="obj">Compare with this instance.</param>
        /// <returns>
        /// <see langword="true"/> if both instances represent the same;
        /// <see langword="false"/> otherwise.
        /// </returns>
        public override bool Equals(object obj)
        {
            return obj is MaterialType type && Equals(type);
        }

        /// <summary>
        /// Get if two <see cref="MaterialType"/> represent the same.
        /// </summary>
        /// <param name="obj">Compare with this instance.</param>
        /// <returns>
        /// <see langword="true"/> if both instances represent the same;
        /// <see langword="false"/> otherwise.
        /// </returns>
        public bool Equals(MaterialType obj)
        {
            return Lighting == obj.Lighting && AlphaMode == obj.AlphaMode;
        }

        /// <summary>
        /// Compute a hash code for the object.
        /// </summary>
        /// <returns>A hash code for the current object.</returns>
        /// <remarks>
        /// * You should not assume that equal hash codes imply object equality.
        /// * You should never persist or use a hash code outside the application domain in which it was created,
        ///   because the same object may hash differently across application domains, processes, and platforms.
        /// </remarks>
        public override int GetHashCode()
        {
            return (Lighting.GetHashCode()) ^ (397 * AlphaMode.GetHashCode());
        }
    }
}
