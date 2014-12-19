/*
    This file is part of Mitsuba, a physically based rendering system.

    Copyright (c) 2007-2014 by Wenzel Jakob and others.

    Mitsuba is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License Version 3
    as published by the Free Software Foundation.

    Mitsuba is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program. If not, see <http://www.gnu.org/licenses/>.
*/

#include <mitsuba/render/volume.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

#include <mitsuba/spgrid/SPGrid_Allocator.h>
#include <mitsuba/spgrid/SPGrid_Set.h>
#include <mitsuba/spgrid/Load_Helper.h>

// Uncomment to enable nearest-neighbor direction interpolation
//#define VINTERP_NEAREST_NEIGHBOR

// Number of power iteration steps used to find the dominant direction
#define POWER_ITERATION_STEPS 5

using namespace SPGrid;

MTS_NAMESPACE_BEGIN

//TODO: change this stuff to match actual SPGRID
/*!\plugin{spgridvolume}{Sparse Paged Grid-based volume data source}
 * \parameters{
 *     \parameter{filename}{\String}{
 *       Specifies the filename of the volume data file to be loaded
 *     }
 *     \parameter{sendData}{\Boolean}{
 *       When this parameter is set to \code{true}, the implementation will
 *       send all volume data to other network render nodes. Otherwise, they
 *       are expected to have access to an identical volume data file that can be
 *       mapped into memory. \default{\code{false}}
 *     }
 *     \parameter{toWorld}{\Transform}{
 *         Optional linear transformation that should be applied to the data
 *     }
 *     \parameter{min, max}{\Point}{
 *         Optional parameter that can be used to re-scale the data so that
 *         it lies in the bounding box between \code{min} and \code{max}.
 *     }
 * }
 *
 * This class implements access to memory-mapped volume data stored on a
 * 3D grid using a simple binary exchange format.
 * The format uses a little endian encoding and is specified as
 * follows:\vspace{3mm}
 *
 * \begin{center}
 * \begin{tabular}{>{\bfseries}p{2cm}p{11cm}}
 * \toprule
 * Position & Content\\
 * \midrule
 * Bytes 1-3&   ASCII Bytes '\code{V}', '\code{O}', and '\code{L}' \\
 * Byte  4&     File format version number (currently 3)\\
 * Bytes 5-8&   Encoding identifier (32-bit integer). The following
 * choices are available:
 * \begin{enumerate}[1.]
 * \item Dense \code{float32}-based representation
 * \item Dense \code{float16}-based representation (\emph{currently not supported by this implementation})
 * \item Dense \code{uint8}-based representation (The range 0..255 will be mapped to 0..1)
 * \item Dense quantized directions. The directions are stored in spherical
 * coordinates with a total storage cost of 16 bit per entry.
 * \end{enumerate}\\
 * Bytes 9-12 &  Number of cells along the X axis (32 bit integer)\\
 * Bytes 13-16 &  Number of cells along the Y axis (32 bit integer)\\
 * Bytes 17-20 &  Number of cells along the Z axis (32 bit integer)\\
 * Bytes 21-24 &  Number of channels (32 bit integer, supported values: 1 or 3)\\
 * Bytes 25-48 &  Axis-aligned bounding box of the data stored in single
 *                precision (order: xmin, ymin, zmin, xmax, ymax, zmax)\\
 * Bytes 49-*  &  Binary data of the volume stored in the specified encoding.
 *                The data are ordered so that the following C-style indexing
 *                operation makes sense after the file has been mapped into memory:\newline
 *                   \ \ \code{data[((zpos*yres + ypos)*xres + xpos)*channels + chan]}\newline
 *                where \code{(xpos, ypos, zpos, chan)} denotes the lookup location.\\
 *
 * \bottomrule
 * \end{tabular}
 * \end{center}
 *
 * Note that Mitsuba expects that entries in direction volumes are either
 * zero or valid unit vectors.
 *
 * When using this data source to represent floating point density volumes,
 * please ensure that the values are all normalized to lie in the
 * range $[0, 1]$---otherwise, the Woodcock-Tracking integration method in
 * \pluginref{heterogeneous} will produce incorrect results.
 */
class SPGridDataSource : public VolumeDataSource {
public:

    typedef struct Foo_struct {
        float x, y, z;
        unsigned flags;
    } Foo;
    typedef SPGrid_Allocator<Foo,3> Foo_Allocator;
    typedef Load_Helper<Foo,3> Loader;
    typedef SPGrid_Allocator<Foo,3>::Array<>::mask Foo_Mask;
    typedef SPGrid_Allocator<Foo,3>::Array<float>::type Data_array_type;
    typedef SPGrid_Allocator<Foo,3>::Array<const float>::type Const_data_array_type;
    typedef SPGrid_Allocator<Foo,3>::Array<unsigned>::type Flags_array_type;
    typedef SPGrid_Set<Flags_array_type> Flags_set_type;
    typedef std_array<int,3> Vec3i;
    typedef std_array<float,3> Vec3f;

	SPGridDataSource(const Properties &props)
		: VolumeDataSource(props),
          alloc(Loader::Load_Allocator(props.getString("blockfile"))),
          d1(alloc.Get_Array(&Foo::x)),
          flags(alloc.Get_Array(&Foo::flags)),
          flag_set(flags)
    {
		m_volumeToWorld = props.getTransform("toWorld", Transform());
		if (props.hasProperty("min") && props.hasProperty("max")) {
			/* Optionally allow to use an AABB other than
			   the one specified by the grid file */
			m_dataAABB.min = props.getPoint("min");
			m_dataAABB.max = props.getPoint("max");
		} else {
            std_array<unsigned, 3> size = alloc.Size();
            //TODO: don't hardcode these
            m_dataAABB.min = Point(-0.5, -0.5, -0.195312);
            m_dataAABB.max = Point( 0.5,  0.5,  0.195312);
        }
        m_res[0] = alloc.xsize;
        m_res[1] = alloc.ysize;
        m_res[2] = alloc.zsize;

		/**
		 * When 'sendData' is set to false, only the filename
		 * is transmitted. A following unserialization of the
		 * stream causes the implementation to then look for
		 * the file (which had better exist if unserialization
		 * occurs on a remote machine).
		 */
        m_sendData = props.getBoolean("sendData", false);

		Loader::Load_Mask(alloc, flag_set, props.getString("blockfile"));
        Loader::Load_Data(alloc, &Foo::flags, flag_set, props.getString("flagfile")); //now it's some type issue?
        Loader::Load_Data(alloc, &Foo::x, flag_set, props.getString("densityfile"));
    }

    //TODO: probably won't support this? or at least force local load?
	SPGridDataSource(Stream *stream, InstanceManager *manager)
			: VolumeDataSource(stream, manager),
              alloc(Loader::Load_Allocator("GIANT FAILURE!")),
              d1(alloc.Get_Array(&Foo::x)),
              flags(alloc.Get_Array(&Foo::flags)),
              flag_set(flags)
    {
        /*
		m_volumeToWorld = Transform(stream);
		m_dataAABB = AABB(stream);
		m_sendData = stream->readBool();
		if (m_sendData) {
			m_volumeType = (EVolumeType) stream->readInt();
			m_res = Vector3i(stream);
			m_channels = stream->readInt();
			m_filename = stream->readString();
			size_t volumeSize = getVolumeSize();
			m_data = new uint8_t[volumeSize];
			stream->read(m_data, volumeSize);
		} else {
			fs::path filename = stream->readString();
			loadFromFile(filename);
		}
		configure();
        */
	}

	virtual ~SPGridDataSource() {
        printf("\nSummary:\n  min: (%4d, %4d, %4d)\n  max: (%4d, %4d, %4d)\n",
            minX, minY, minZ,
            maxX, maxY, maxZ);
        printf(" Caused By:\n");
        printf("  (%f, %f, %f)\n", nix.x, nix.y, nix.z);
        printf("  (%f, %f, %f)\n", niy.x, niy.y, niy.z);
        printf("  (%f, %f, %f)\n", niz.x, niz.y, niz.z);
        printf("  (%f, %f, %f)\n", xix.x, xix.y, xix.z);
        printf("  (%f, %f, %f)\n", xiy.x, xiy.y, xiy.z);
        printf("  (%f, %f, %f)\n", xiz.x, xiz.y, xiz.z);
	}

	size_t getVolumeSize() const {
        return size_t(alloc.Padded_Volume());
	}
    
    //TODO: we probably won't support this??? idk
	void serialize(Stream *stream, InstanceManager *manager) const {
		VolumeDataSource::serialize(stream, manager);
        /*
		m_volumeToWorld.serialize(stream);
		m_dataAABB.serialize(stream);
		stream->writeBool(m_sendData);

		if (m_sendData) {
			stream->writeInt(m_volumeType);
			m_res.serialize(stream);
			stream->writeInt(m_channels);
			stream->writeString(m_filename.string());
			stream->write(m_data, getVolumeSize());
		} else {
			stream->writeString(m_filename.string());
		}
        */
	}

    void loadFromFile(const fs::path &filepath) {
        printf("Cannot load from file\n");
        //TODO
    }

    //TODO: fix me for SPGrid
	void configure() {
		Vector extents(m_dataAABB.getExtents());
		m_worldToVolume = m_volumeToWorld.inverse();
		m_worldToGrid = Transform::scale(Vector(
				(m_res[0] - 1) / extents[0],
				(m_res[1] - 1) / extents[1],
				(m_res[2] - 1) / extents[2])
			) * Transform::translate(-Vector(m_dataAABB.min)) * m_worldToVolume;
		m_stepSize = std::numeric_limits<Float>::infinity();
		for (int i=0; i<3; ++i)
			m_stepSize = 0.5f * std::min(m_stepSize, extents[i] / (Float) (m_res[i]-1));
		m_aabb.reset();
		for (int i=0; i<8; ++i)
			m_aabb.expandBy(m_volumeToWorld(m_dataAABB.getCorner(i)));

        /* Precompute cosine and sine lookup tables */
		/*
        for (int i=0; i<255; i++) {
			Float angle = (float) i * ((float) M_PI / 255.0f);
			m_cosPhi[i] = std::cos(2.0f * angle);
			m_sinPhi[i] = std::sin(2.0f * angle);
			m_cosTheta[i] = std::cos(angle);
			m_sinTheta[i] = std::sin(angle);
			m_densityMap[i] = i/255.0f;
		}
		m_cosPhi[255] = m_sinPhi[255] = 0;
		m_cosTheta[255] = m_sinTheta[255] = 0;
		m_densityMap[255] = 1.0f;
        */

        minX = minY = minZ = 0x7FFFFFFF;
        maxX = maxY = maxZ = 0x80000000;
	}

	/**
	 * This is needed since Mitsuba might be
	 * compiled with either single/double precision
	 */
	struct float3 {
		float value[3];

		inline float3() { }

		inline float3(float a, float b, float c) {
			value[0] = a; value[1] = b; value[2] = c;
		}

		inline explicit float3(double a, double b, double c) {
			value[0] = (float) a; value[1] = (float) b; value[2] = (float) c;
		}

		inline float3 operator*(Float v) const {
			return float3((float) (value[0]*v), (float) (value[1]*v), (float) (value[2]*v));
		}

		inline float3 operator+(const float3 &f2) const {
			return float3(value[0]+f2.value[0], value[1]+f2.value[1], value[2]+f2.value[2]);
		}

		inline Spectrum toSpectrum() const {
			Spectrum result;
			result.fromLinearRGB(value[0], value[1], value[2]);
			return result;
		}

		inline Vector toVector() const {
			return Vector(value[0], value[1], value[2]);
		}

		float operator[](int i) const {
			return value[i];
		}

		inline Matrix3x3 tensor() const {
			return Matrix3x3(
				value[0]*value[0], value[0]*value[1], value[0]*value[2],
				value[1]*value[0], value[1]*value[1], value[1]*value[2],
				value[2]*value[0], value[2]*value[1], value[2]*value[2]
			);
		}
	};

    //TODO: replace me with SPGrid accessing stuff!
    //      it looks like d000-d111 are the grid at the point and the 8 surrounding grids?
    //      need to determine dxxx map to which directions
    //      or just grab the closest thing and call it good? idk
	int minX, minY, minZ, maxX, maxY, maxZ;
    Point nix, niy, niz, xix, xiy, xiz;
    Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);

		//This is just the basic implementation - grab nearest neighbor?
		int x,y,z;
		x = math::floorToInt(p.x);
		y = math::floorToInt(p.y);
		z = math::floorToInt(p.z);
        std_array<int,3> coord(x, y, z);
        if (flag_set.Is_Set(coord, 0xFFFFFFFFU)) {
		    return d1(coord);
        }
        return 0.0f;

		//The following stuff does some nice interpolation but we can look at that later
		/*
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return 0;

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
				const float *floatData = (float *) m_data;
				const Float
					d000 = floatData[(z1*m_res.y + y1)*m_res.x + x1],
					d001 = floatData[(z1*m_res.y + y1)*m_res.x + x2],
					d010 = floatData[(z1*m_res.y + y2)*m_res.x + x1],
					d011 = floatData[(z1*m_res.y + y2)*m_res.x + x2],
					d100 = floatData[(z2*m_res.y + y1)*m_res.x + x1],
					d101 = floatData[(z2*m_res.y + y1)*m_res.x + x2],
					d110 = floatData[(z2*m_res.y + y2)*m_res.x + x1],
					d111 = floatData[(z2*m_res.y + y2)*m_res.x + x2];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			case EUInt8: {
				const Float
					d000 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x1]],
					d001 = m_densityMap[m_data[(z1*m_res.y + y1)*m_res.x + x2]],
					d010 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x1]],
					d011 = m_densityMap[m_data[(z1*m_res.y + y2)*m_res.x + x2]],
					d100 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x1]],
					d101 = m_densityMap[m_data[(z2*m_res.y + y1)*m_res.x + x2]],
					d110 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x1]],
					d111 = m_densityMap[m_data[(z2*m_res.y + y2)*m_res.x + x2]];

				return ((d000*_fx + d001*fx)*_fy +
						(d010*_fx + d011*fx)*fy)*_fz +
					   ((d100*_fx + d101*fx)*_fy +
						(d110*_fx + d111*fx)*fy)*fz;
			}
			default:
				return 0.0f;
		}
		*/
	}

    //TODO: are we supporting this? probably just force return Spectrum(0.0f) for now...
	Spectrum lookupSpectrum(const Point &_p) const {
		return Spectrum(0.0f);
		/*
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return Spectrum(0.0f);

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1,
				_fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

		switch (m_volumeType) {
			case EFloat32: {
				const float3 *spectrumData = (float3 *) m_data;
				const float3
					&d000 = spectrumData[(z1*m_res.y + y1)*m_res.x + x1],
					&d001 = spectrumData[(z1*m_res.y + y1)*m_res.x + x2],
					&d010 = spectrumData[(z1*m_res.y + y2)*m_res.x + x1],
					&d011 = spectrumData[(z1*m_res.y + y2)*m_res.x + x2],
					&d100 = spectrumData[(z2*m_res.y + y1)*m_res.x + x1],
					&d101 = spectrumData[(z2*m_res.y + y1)*m_res.x + x2],
					&d110 = spectrumData[(z2*m_res.y + y2)*m_res.x + x1],
					&d111 = spectrumData[(z2*m_res.y + y2)*m_res.x + x2];

				return (((d000*_fx + d001*fx)*_fy +
						 (d010*_fx + d011*fx)*fy)*_fz +
						((d100*_fx + d101*fx)*_fy +
						 (d110*_fx + d111*fx)*fy)*fz).toSpectrum();
				}
			case EUInt8: {
				const float3
					d000 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x1)+2]]),
					d001 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y1)*m_res.x + x2)+2]]),
					d010 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x1)+2]]),
					d011 = float3(
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z1*m_res.y + y2)*m_res.x + x2)+2]]),
					d100 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x1)+2]]),
					d101 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y1)*m_res.x + x2)+2]]),
					d110 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x1)+2]]),
					d111 = float3(
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+0]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+1]],
						m_densityMap[m_data[3*((z2*m_res.y + y2)*m_res.x + x2)+2]]);

				return (((d000*_fx + d001*fx)*_fy +
						 (d010*_fx + d011*fx)*fy)*_fz +
						((d100*_fx + d101*fx)*_fy +
						 (d110*_fx + d111*fx)*fy)*fz).toSpectrum();

				}
			default: return Spectrum(0.0f);
		}
		*/
	}

    //TODO: same sort of thing as with lookupSpectrum - could be done but not necessary for what we're doing??
	Vector lookupVector(const Point &_p) const {
		return Vector(0.0f);
		/*
		const Point p = m_worldToGrid.transformAffine(_p);
		const int x1 = math::floorToInt(p.x),
			  y1 = math::floorToInt(p.y),
			  z1 = math::floorToInt(p.z),
			  x2 = x1+1, y2 = y1+1, z2 = z1+1;

		if (x1 < 0 || y1 < 0 || z1 < 0 || x2 >= m_res.x ||
		    y2 >= m_res.y || z2 >= m_res.z)
			return Vector(0.0f);

		const Float fx = p.x - x1, fy = p.y - y1, fz = p.z - z1;
		Vector value;
		*/
		//#if defined(VINTERP_NEAREST_NEIGHBOR)
			/* Nearest neighbor */
		/*
			switch (m_volumeType) {
				case EFloat32: {
					const float3 *vectorData = (float3 *) m_data;
					value = vectorData[
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2)].toVector();
					}
					break;
				case EQuantizedDirections: {
					value = lookupQuantizedDirection(
						(((fz < .5) ? z1 : z2) * m_res.y +
						((fy < .5) ? y1 : y2)) * m_res.x +
						((fx < .5) ? x1 : x2));
					}
					break;
				default:
					return Vector(0.0f);
			}
		#else
			Float _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

			Matrix3x3 tensor(0.0f);
			switch (m_volumeType) {
				case EFloat32: {
						const float3 *vectorData = (float3 *) m_data;
						for (int k=0; k<8; ++k) {
							uint32_t index = (((k & 4) ? z2 : z1) * m_res.y +
								((k & 2) ? y2 : y1)) * m_res.x + ((k & 1) ? x2 : x1);
							Float factor = ((k & 1) ? fx : _fx) * ((k & 2) ? fy : _fy)
								* ((k & 4) ? fz : _fz);
							Vector d = vectorData[index].toVector();
							tensor(0, 0) += factor * d.x * d.x;
							tensor(0, 1) += factor * d.x * d.y;
							tensor(0, 2) += factor * d.x * d.z;
							tensor(1, 1) += factor * d.y * d.y;
							tensor(1, 2) += factor * d.y * d.z;
							tensor(2, 2) += factor * d.z * d.z;
						}
					}
					break;
				case EQuantizedDirections: {
						for (int k=0; k<8; ++k) {
							uint32_t index = (((k & 4) ? z2 : z1) * m_res.y +
								((k & 2) ? y2 : y1)) * m_res.x + ((k & 1) ? x2 : x1);
							Float factor = ((k & 1) ? fx : _fx) * ((k & 2) ? fy : _fy)
								* ((k & 4) ? fz : _fz);
							Vector d = lookupQuantizedDirection(index);
							tensor(0, 0) += factor * d.x * d.x;
							tensor(0, 1) += factor * d.x * d.y;
							tensor(0, 2) += factor * d.x * d.z;
							tensor(1, 1) += factor * d.y * d.y;
							tensor(1, 2) += factor * d.y * d.z;
							tensor(2, 2) += factor * d.z * d.z;
						}
					}
					break;
				default:
					return Vector(0.0f);
			}

			tensor(1, 0) = tensor(0, 1);
			tensor(2, 0) = tensor(0, 2);
			tensor(2, 1) = tensor(1, 2);

			if (tensor.isZero())
				return Vector(0.0f);

#if 0
			Float lambda[3];
			eig3_noniter(tensor, lambda);
			value = tensor.col(0);
			Float specularity = 1-lambda[1]/lambda[0];

#else
			* Square the structure tensor for faster convergence */
			/*
			tensor *= tensor;

			const Float invSqrt3 = 0.577350269189626f;
			value = Vector(invSqrt3, invSqrt3, invSqrt3);

			* Determine the dominant eigenvector using
			   a few power iterations */
			/*
			for (int i=0; i<POWER_ITERATION_STEPS-1; ++i)
				value = normalize(tensor * value);
			value = tensor * value;
#endif

		#endif

		if (!value.isZero())
			return normalize(m_volumeToWorld(value));
		else
			return Vector(0.0f);
		*/
	}

    //TODO: adjust these to show what we actually support i.e. true/false/false if we don't implement Spectrum/Vector stuff
	bool supportsFloatLookups() const {return true;}
	bool supportsSpectrumLookups() const { return false; }
	bool supportsVectorLookups() const { return false; }
	/*	
	bool supportsFloatLookups() const { return m_channels == 1; }
	bool supportsSpectrumLookups() const { return m_channels == 3; }
	bool supportsVectorLookups() const { return m_channels == 3; }
	*/	
	Float getStepSize() const { return m_stepSize; }

	Float getMaximumFloatValue() const {
		return 1.0f;
	}

	std::string toString() const {
		std::ostringstream oss;
		oss << "SPGridVolume[" << endl
			<< "  res = " << m_res.toString() << "," << endl
			<< "  aabb = " << m_dataAABB.toString() << endl
			<< "]";
		return oss.str();
	}

	MTS_DECLARE_CLASS()
protected:
    /*
	FINLINE Vector lookupQuantizedDirection(size_t index) const {
		uint8_t theta = m_data[2*index], phi = m_data[2*index+1];
		return Vector(
			m_cosPhi[phi] * m_sinTheta[theta],
			m_sinPhi[phi] * m_sinTheta[theta],
			m_cosTheta[theta]
		);
	}
    */

protected:
    Foo_Allocator alloc;
    Data_array_type d1;
    Flags_array_type flags;
    Flags_set_type flag_set;
	bool m_sendData;
	Vector3i m_res;
	Transform m_worldToGrid;
	Transform m_worldToVolume;
	Transform m_volumeToWorld;
	Float m_stepSize;
	AABB m_dataAABB;
	//Float m_cosTheta[256], m_sinTheta[256];
	//Float m_cosPhi[256], m_sinPhi[256];
	//Float m_densityMap[256];
};

MTS_IMPLEMENT_CLASS_S(SPGridDataSource, false, VolumeDataSource);
MTS_EXPORT_PLUGIN(SPGridDataSource, "Sparse Paged Grid data source");
MTS_NAMESPACE_END