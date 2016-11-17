#pragma once 
#include <memory>
#include <CPT/engine/data_pool/pipeline_schema.hpp>
#include <CPT/engine/data_pool/data_paths_pool.hpp>
#include <CPT/engine/data_pool/shared_object_manager.hpp>
#include <CPT/engine/data_pool/monitor.hpp>
#include <CPT/engine/data_pool/parallel_thread_pool.hpp>
#include <CCD/utility/language.hpp>
#include <AGO/engine/data_pool/fastq.hpp>
#include <AGO/engine/data_pool/genome.hpp>
#include <AGO/engine/data_pool/tailor.hpp>

namespace ago {
namespace engine {

auto& monitor_output_ = cpt::msg; //cpt::msg;
std::ofstream mointor_outfile_;

class DataPoolImpl 
// [Core]
: public cpt::engine::data_pool::PipelineSchema // input
, public cpt::engine::data_pool::DataPathsPool
, public cpt::engine::data_pool::ParallelThreadPool
, public cpt::engine::data_pool::SharedObjectManager
, public ago::engine::data_pool::FastqImpl
, public ago::engine::data_pool::GenomeImpl
, public ago::engine::data_pool::Tailor
// [Opts]
{
    cpt::engine::data_pool::Monitor< decltype( monitor_output_ )> monitor_;

  public:

    template<class OP>
    DataPoolImpl(OP& op)
    : cpt::engine::data_pool::PipelineSchema( op.pipeline_schema_stream_ )
    , cpt::engine::data_pool::ParallelThreadPool( 20 )
    , FastqImpl( *this )
    , GenomeImpl( *this )
    , Tailor( *this )
    , monitor_( monitor_output_, mointor_outfile_ )
    {
        monitor_.set_output_dir( this->output_dir().string() );
    }

    auto& monitor()
    {
        return monitor_;
    }

    DEFAULT_MOVE(DataPoolImpl);
    DISABLE_COPY(DataPoolImpl);
};

using DataPool = DataPoolImpl;
using DataPoolPtr = std::unique_ptr<DataPool>;

} //engine
} //ago
