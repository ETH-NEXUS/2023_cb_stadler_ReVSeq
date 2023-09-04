import {TableDataSampleCounts, TableDataMetadata} from 'src/models/tables'
import {SampleCount, Metadata} from 'src/models/core'

export const createSampleCountsData = (sampleCounts: SampleCount[]): TableDataSampleCounts[] => {
  return sampleCounts.map(sampleCount => {
    return {
      pseudoanonymized_id: sampleCount.sample
        ? sampleCount.sample.pseudoanonymized_id
        : sampleCount.pseudoanonymized_id,
      substrain: sampleCount.substrain.name,
      strain: sampleCount.substrain.strain.name,
      aligned: sampleCount.aligned,
      length: sampleCount.length,
      rpkm: sampleCount.rpkm,
      rpkm_proportions: sampleCount.rpkm_proportions,
      //  normcounts: sampleCount.normcounts,
      outlier: sampleCount.outlier,
    }
  })
}

export const createMetadataRows = (metadata: Metadata[]): TableDataMetadata[] => {
  return metadata.map(m => {
    return {
      pseudoanonymized_id: m.sample.pseudoanonymized_id,
      prescriber: m.prescriber,
      well: m.well.location,
      order_date: m.order_date,
      ent_date: m.ent_date,
      treatment_type: m.treatment_type,
    }
  })
}
