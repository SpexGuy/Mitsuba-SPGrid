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

//TODO: remove this line (forces debug compilation)
#undef MTS_NDEBUG

#include <mitsuba/render/volume.h>
#include <mitsuba/core/mstream.h>
#include <mitsuba/core/fresolver.h>
#include <mitsuba/core/properties.h>

#include <mitsuba/spgrid/SPGrid_Allocator.h>
#include <mitsuba/spgrid/SPGrid_Set.h>
#include <mitsuba/spgrid/Load_Helper.h>

// Uncomment to enable nearest-neighbor direction interpolation
//#define VINTERP_NEAREST_NEIGHBOR
#define VINTERP_LOOKUP_FLOAT
//#define VINTERP_LOOKUP_SPECTRUM

// Number of power iteration steps used to find the dominant direction
//#define POWER_ITERATION_STEPS 5

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

    struct SPGrid {
    public:
        Foo_Allocator alloc;
        Data_array_type d0;
        Data_array_type d1;
        Data_array_type d2;
        Flags_array_type flags;
        Flags_set_type flag_set;
        unsigned entry_size;

        SPGrid(Foo_Allocator &&alloc_in, unsigned entry_size_in)
            : alloc(std::move(alloc_in)),
              d0(alloc.Get_Array(&Foo::x)),
              d1(alloc.Get_Array(&Foo::y)),
              d2(alloc.Get_Array(&Foo::z)),
              flags(alloc.Get_Array(&Foo::flags)),
              flag_set(flags)
        {
            entry_size = entry_size_in;
        }
    };

	SPGridDataSource(const Properties &props)
        : VolumeDataSource(props)
    {
        bool error = false;
        
        std::vector<std::string> blockfiles;
        std::vector<std::string> flagfiles;
        std::vector<std::string> channel0files;
        std::vector<std::string> channel1files;
        std::vector<std::string> channel2files;

        std::string blocks = props.getString("blockfile");
        std::string flags = props.getString("flagfile");
        std::string channel0 = props.getString("channel0file");

        parseFileList(blockfiles, "blockfile", blocks);
        parseFileList(flagfiles, "flagfile", flags);
        parseFileList(channel0files, "channel0file", channel0);

        m_numChannels = 1;
        if (props.hasProperty("channel1file") && props.hasProperty("channel2file")) {
            m_numChannels = 3;
            std::string channel1 = props.getString("channel1file");
            std::string channel2 = props.getString("channel2file");
            parseFileList(channel1files, "channel1file", channel1);
            parseFileList(channel2files, "channel2file", channel2);
        } else if (props.hasProperty("channel1file")) {
            Log(EError, "channel1file specified but not channel2file!");
            error = true;
        } else if (props.hasProperty("channel2file")) {
            Log(EError, "channel2file specified but not channel1file!");
            error = true;
        }

        m_nGrids = blockfiles.size();
        if (m_nGrids == 0) {
            Log(EError, "No grids specified! Use the blockfile, flagfile, and channel0file attributes.");
            error = true;
        } else {
            if (flagfiles.size() != m_nGrids) {
                Log(EError, "Wrong number of flagfiles! (%d blockfiles, %d flagfiles)", m_nGrids, flagfiles.size());
                error = true;
            }
            if (channel0files.size() != m_nGrids) {
                Log(EError, "Wrong number of channel0files! (%d blockfiles, %d channel0files)", m_nGrids, channel0files.size());
                error = true;
            }
            if (m_numChannels == 3) {
                if (channel1files.size() != m_nGrids) {
                    Log(EError, "Wrong number of channel1files! (%d blockfiles, %d channel1files)", m_nGrids, channel1files.size());
                    error = true;
                }
                if (channel2files.size() != m_nGrids) {
                    Log(EError, "Wrong number of channel2files! (%d blockfiles, %d channel2files)", m_nGrids, channel2files.size());
                    error = true;
                }
            }
        }

        if (error) {
            Log(EError, "At least one thing went wrong loading the SPGrid. Setting number of grids to 0.");
            m_nGrids = 0;
        }

        if (m_nGrids) {
            m_grids = new SPGrid*[m_nGrids];
            for (unsigned c = 0; c < m_nGrids; c++) {
                m_grids[c] = new SPGrid(std::move(Loader::Load_Allocator(blockfiles[c])), c);
                Loader::Load_Mask(m_grids[c]->alloc, m_grids[c]->flag_set, blockfiles[c]);
                Loader::Load_Data(m_grids[c]->alloc, &Foo::flags, m_grids[c]->flag_set, flagfiles[c]);
                Loader::Load_Data(m_grids[c]->alloc, &Foo::x, m_grids[c]->flag_set, channel0files[c]);
                if (m_numChannels == 3) {
                    Loader::Load_Data(m_grids[c]->alloc, &Foo::y, m_grids[c]->flag_set, channel1files[c]);
                    Loader::Load_Data(m_grids[c]->alloc, &Foo::z, m_grids[c]->flag_set, channel2files[c]);
                }
            }
        }

        m_volumeToWorld = props.getTransform("toWorld", Transform());

        if (props.hasProperty("min") && props.hasProperty("max")) {
			/* Optionally allow to use an AABB other than
			   the one specified by the grid file */
			m_dataAABB.min = props.getPoint("min");
			m_dataAABB.max = props.getPoint("max");
		} else {
            //TODO: don't hardcode these
            m_dataAABB.min = Point(-0.5, -0.5, -0.195312);
            m_dataAABB.max = Point( 0.5,  0.5,  0.195312);
        }

        if (m_nGrids) {
            m_res[0] = m_grids[0]->alloc.xsize;
            m_res[1] = m_grids[0]->alloc.ysize;
            m_res[2] = m_grids[0]->alloc.zsize;
        } else {
            m_res[0] = 1;
            m_res[1] = 1;
            m_res[2] = 1;
        }

        printf("Loaded %d SPGrids\n", m_nGrids);
        
        //TODO: validate grid sizes

		/**
		 * When 'sendData' is set to false, only the filename
		 * is transmitted. A following unserialization of the
		 * stream causes the implementation to then look for
		 * the file (which had better exist if unserialization
		 * occurs on a remote machine).
		 */
        m_sendData = props.getBoolean("sendData", false);
        
    }

    //TODO: probably won't support this? or at least force local load?
	SPGridDataSource(Stream *stream, InstanceManager *manager)
        : VolumeDataSource(stream, manager)
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
        for (unsigned c = 0; c < m_nGrids; c++)
            delete m_grids[c];
        delete [] m_grids;
	}

/*
	size_t getVolumeSize() const {
        return size_t(alloc.Padded_Volume());
	}
*/
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

    Float lookupFloat(const Point &_p) const {
		const Point p = m_worldToGrid.transformAffine(_p);

		int x1 = math::floorToInt(p.x),
            y1 = math::floorToInt(p.y),
            z1 = math::floorToInt(p.z);
		
        if (x1 < 0 || y1 < 0 || z1 < 0)
            return 0;
        if (x1 >= m_res.x || y1 >= m_res.y || z1 >= m_res.z)
			return 0;
#ifdef VINTERP_LOOKUP_FLOAT
        float d[8];
        
        const unsigned resolution = getInitialGridValue(x1, y1, z1, d[0]);
        const int x2 = x1 + resolution,
                  y2 = y1 + resolution,
                  z2 = z1 + resolution;

        getGridValue(x2, y1, z1, d[1]);
        getGridValue(x1, y2, z1, d[2]);
        getGridValue(x2, y2, z1, d[3]);
        getGridValue(x1, y1, z2, d[4]);
        getGridValue(x2, y1, z2, d[5]);
        getGridValue(x1, y2, z2, d[6]);
        getGridValue(x2, y2, z2, d[7]);

        const Float fx = (p.x - x1)/resolution, fy = (p.y - y1)/resolution, fz = (p.z - z1)/resolution,
                   _fx = 1.0f - fx, _fy = 1.0f - fy, _fz = 1.0f - fz;

        return ((d[0]*_fx + d[1]*fx)*_fy +
                (d[2]*_fx + d[3]*fx)*fy)*_fz +
               ((d[4]*_fx + d[5]*fx)*_fy +
                (d[6]*_fx + d[7]*fx)*fy)*fz;
#else
        float value;
        getGridValue(x1, y1, z1, value);
        return value;
#endif
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
	bool supportsFloatLookups() const { return true; }
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
    // trim from start
    static inline std::string &ltrim(std::string &s) {
            s.erase(s.begin(), std::find_if(s.begin(), s.end(), std::not1(std::ptr_fun<int, int>(std::isspace))));
            return s;
    }

    // trim from end
    static inline std::string &rtrim(std::string &s) {
            s.erase(std::find_if(s.rbegin(), s.rend(), std::not1(std::ptr_fun<int, int>(std::isspace))).base(), s.end());
            return s;
    }

    // trim from both ends
    static inline std::string &trim(std::string &s) {
            return ltrim(rtrim(s));
    }

    FINLINE void parseFileList(std::vector<std::string> &files, const std::string &name, const std::string &str) {
        std::stringstream ss(str);
        std::string item;
        while (std::getline(ss, item, ',')) {
            item = trim(item);
            if (item.empty())
                Log(EWarn, "Item %s contains a blank file at index %d.", name.c_str(), files.size());
            else
                files.push_back(item);
        }
    }

    //shifting the coordinates by c to account for the resolution difference in the lower grids.
    //TODO: optimize based on number of trailing zeros
    FINLINE unsigned getInitialGridValue(int &x, int &y, int &z, float &val) const {
        for (unsigned c = 0; c < m_nGrids; c++) {
            std_array<int,3> coord((x>>c), (y>>c), (z>>c));
            if (m_grids[c]->flag_set.Is_Set(coord, 0xFFFFFFFFU)) {
                val = m_grids[c]->d0(coord);
                x &= 0xFFFFFFFF << c; // zero the lowest [c] bits in x, y, and z to align it with the grid
                y &= 0xFFFFFFFF << c;
                z &= 0xFFFFFFFF << c;
                return 1 << c;
            }
        }
        val = 0.0f;
        return 1;
    }

    FINLINE void getGridValue(int x, int y, int z, float &val) const {
        for (unsigned c = 0; c < m_nGrids; c++) {
            std_array<int,3> coord((x>>c), (y>>c), (z>>c));
            if (m_grids[c]->flag_set.Is_Set(coord, 0xFFFFFFFFU)) {
                val = m_grids[c]->d0(coord);
                return;
            }
        }
        val = 0.0f;
    }

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
    SPGrid **m_grids;
    unsigned m_nGrids;
    unsigned m_numChannels;
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
