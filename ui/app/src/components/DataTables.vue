<script setup lang="ts">
import {columnsSampleCount, columnsMetadata, SampleCountsTableColumns} from 'components/data'
import {createSampleCountsData, createMetadataRows} from 'components/helpers'
import {useCoreStore} from 'stores/core'
import {useI18n} from 'vue-i18n'
import {computed} from 'vue'

const {t} = useI18n()

const coreStore = useCoreStore()

const aggregate = () => {
  try {
    coreStore.aggregate = !coreStore.aggregate
    coreStore.toggleAggregate()
  } catch (err) {
    console.error(err)
  }
}
const filterColumnNames = computed(() => {
  if (coreStore.aggregate) {
    return columnsSampleCount.filter((c: SampleCountsTableColumns) => c.name !== 'substrain')
  }
  return columnsSampleCount
})
</script>

<template>
  <div class="q-pa-md" v-if="coreStore.sampleCounts.length > 0">
    <h3 class="text-h4 text-center q-pb-md">Counts</h3>
    <div class="tw-flex tw-flex-wrap tw-mb-20 tw-justify-center"></div>

    <div class="tw-mb-4">
      <q-btn
        icon="query_stats"
        :label="t('label.aggregate')"
        color="primary"
        class="tw-mr-2"
        @click="aggregate"></q-btn>
      <q-btn color="primary" icon="science" :label="t('label.raw_data')" @click="aggregate"></q-btn>
    </div>
    <q-table
      :rows-per-page-options="[10, 20, 50, 100]"
      rows-per-page="20"
      flat
      bordered
      :rows="createSampleCountsData(coreStore.tableData)"
      :columns="filterColumnNames"
      row-key="name"></q-table>
  </div>
  <div class="q-pa-md" v-if="coreStore.metadata.length > 0">
    <h3 class="text-h4 text-center q-pb-md">Metadata</h3>
    <q-table
      :rows-per-page-options="[50, 100]"
      rows-per-page="50"
      flat
      bordered
      :rows="createMetadataRows(coreStore.metadata)"
      :columns="columnsMetadata"
      row-key="name"></q-table>
  </div>
</template>

<style scoped>
.point {
  cursor: pointer;
}
.point:hover {
  background-color: darkslateblue;
  cursor: pointer;
}
</style>
